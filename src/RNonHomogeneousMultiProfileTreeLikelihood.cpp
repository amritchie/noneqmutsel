

/** File: RNonHomogeneousMultiProfileTreeLikelihood
 **
 ** Very much in progress, sorry.
 **
 ** This class is intended to give the likelihood of the alignment under a model where each site
 ** evolves under a mixture of tree-wide evolutionary histories. 
 ** Each history should be a SubstitutionModelSet allowing different Mutation Selection Models
 ** to be applied following each breakpoint.
 ** The breakpoints will be the same for all sites.
 ** Most of this stuff is identical to functions from other classes in Bio++; the class hierarchy
 ** seems to have to be rewritten quite a ways back to allow for this kind of branch-site modelling.
 **/

#include <algorithm>
#include <vector>

#include <Bpp/Phyl/Model/MixtureOfSubstitutionModels.h>
#include <Bpp/Phyl/PatternTools.h>
#include "/home/aritchie/src/nemss2/src/RNonHomogeneousMultiProfileTreeLikelihood.h"
//#include "/home/aritchie/src/nemss2/src/NonHomogeneousSiteProfileMixture.h"

using namespace bpp;

using namespace std;

/****************************************************************************************************/


RNonHomogeneousMultiProfileTreeLikelihood::RNonHomogeneousMultiProfileTreeLikelihood
(const Tree& tree,
 const SiteContainer& sites,  
 vector <SubstitutionModelSet*>& profileSet,
 SimpleDiscreteDistribution* profileFreqs) :
  profileSet_(profileSet),
  profileFreqs_(profileFreqs),
  brLenParameters_(),
  pxy_(),
  dpxy_(),
  d2pxy_(),
  rootFreqs_(),
  nodes_(),
  idToNode_(),
  nbSites_(),
  nbDistinctSites_(),
  nbProfiles_(),
  nbStates_(),
  nbNodes_(),
  minimumBrLen_(),
  maximumBrLen_(),
  brLenConstraint_(),
  root1_(),
  root2_(),
  likelihoodData_(0),
  minusLogLik_(-1.)
{

  AbstractTreeLikelihood::enableDerivatives(true);

  for (auto p : profileSet_) {
    if (!p->isFullySetUpFor(tree))
      throw
	Exception("RNonHomogeneousMultiProfileTreeLikelihood(constructor). Model set is not fully specified.");
  }
    
init_(tree, profileSet_, profileFreqs_);
setData(sites);
}

RNonHomogeneousMultiProfileTreeLikelihood::RNonHomogeneousMultiProfileTreeLikelihood (const RNonHomogeneousMultiProfileTreeLikelihood& lik) :
  AbstractTreeLikelihood(lik),
  profileSet_(lik.profileSet_),
  profileFreqs_(lik.profileFreqs_),
  brLenParameters_(lik.brLenParameters_),
  pxy_(lik.pxy_),
  dpxy_(lik.dpxy_),
  d2pxy_(lik.d2pxy_),
  rootFreqs_(lik.rootFreqs_),
  nodes_(),
  idToNode_(),
  nbSites_(lik.nbSites_),
  nbDistinctSites_(lik.nbDistinctSites_),
  nbProfiles_(lik.nbProfiles_),
  nbStates_(lik.nbStates_),
  nbNodes_(lik.nbNodes_),
  minimumBrLen_(lik.minimumBrLen_),
  maximumBrLen_(lik.maximumBrLen_),
  brLenConstraint_(dynamic_cast<Constraint*>(lik.brLenConstraint_->clone())),
  root1_(lik.root1_),
  root2_(lik.root2_),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  AbstractTreeLikelihood::enableDerivatives(true);
  nodes_ = tree_->getNodes();
  nodes_.pop_back();
  for (unsigned int i = 0; i < nodes_.size(); i++) {
    const Node* node = nodes_[i];
    idToNode_[node->getId()] = node;
  }

  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
}

RNonHomogeneousMultiProfileTreeLikelihood& RNonHomogeneousMultiProfileTreeLikelihood::operator=(const RNonHomogeneousMultiProfileTreeLikelihood& lik){

  AbstractTreeLikelihood::operator=(lik);
  
  profileSet_        = lik.profileSet_;
  profileFreqs_      = lik.profileFreqs_; 
  brLenParameters_   = lik.brLenParameters_;
  pxy_               = lik.pxy_;
  dpxy_              = lik.dpxy_;
  d2pxy_             = lik.d2pxy_;
  rootFreqs_         = lik.rootFreqs_;
  nodes_             = tree_->getNodes();
  nodes_.pop_back(); //Remove the root node (the last added!).  
  nbSites_           = lik.nbSites_;
  nbDistinctSites_   = lik.nbDistinctSites_;
  nbProfiles_         = lik.nbProfiles_;
  nbStates_          = lik.nbStates_;
  nbNodes_           = lik.nbNodes_;
  minimumBrLen_      = lik.minimumBrLen_;
  maximumBrLen_      = lik.maximumBrLen_;
  if (brLenConstraint_.get()) brLenConstraint_.release();
  brLenConstraint_.reset(lik.brLenConstraint_->clone());
  root1_             = lik.root1_;
  root2_             = lik.root2_;

  for (unsigned int i = 0; i < nodes_.size(); i++) {

    const Node * node = nodes_[i];
    idToNode_[node->getId()] = node;

  }
  

  if (likelihoodData_) delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;

}

void RNonHomogeneousMultiProfileTreeLikelihood::init_(const Tree& tree,
						      vector<SubstitutionModelSet*>& modelSet,
						      SimpleDiscreteDistribution* rDist)
{
  TreeTools::checkIds(tree,true);
  tree_ = new TreeTemplate<Node>(tree);
  root1_ = tree_->getRootNode()->getSon(0)->getId();
  root2_ = tree_->getRootNode()->getSon(1)->getId();
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); //Remove the root node
  nbNodes_ = nodes_.size();
  //Build nodes index:
  for (unsigned int i = 0; i < nodes_.size(); i++)
    {
      const Node * node = nodes_[i];
      idToNode_[node->getId()] = node;
    }
  nbProfiles_ = modelSet.size();

  minimumBrLen_ = 0.000001;
  maximumBrLen_ = 10000;
  brLenConstraint_.reset(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true));
  setProfiles(modelSet);

  likelihoodData_ = new DRASRTreeLikelihoodData(tree_,
						nbProfiles_,
						true);
}

void RNonHomogeneousMultiProfileTreeLikelihood::setProfiles(vector <SubstitutionModelSet*>& modelSet) {

  if (data_) {
    //Check Alphabet
    for (auto p :modelSet) {
      if (p->getAlphabet()->getAlphabetType() != data_->getAlphabet()->getAlphabetType())
	throw
	  Exception("RMultiProfileNonHomogeneousTreeLikelihood::setProfiles(). Model alphabet does not match existing data.");
    }
  }

  profileSet_ = modelSet;

  nbStates_ = profileSet_[0]->getNumberOfStates();
  nbProfiles_ = profileSet_.size();

  rootFreqs_.resize(nbProfiles_);

  //Allocate transition probabilities arrays:
  for (unsigned int l = 0; l < nbNodes_; l++) {
    //For each son node, i
    Node * son = nodes_[l];

    VVVdouble* pxy__son = & pxy_[son->getId()];
    pxy__son->resize(nbProfiles_);
    for (unsigned int p = 0; p < nbProfiles_; p++) {
      VVdouble * pxy__son_p = & (*pxy__son)[p];
      pxy__son_p->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++) {
	(*pxy__son_p)[x].resize(nbStates_);
      }
    }
    
    VVVdouble* dpxy__son = & dpxy_[son->getId()];
    dpxy__son->resize(nbProfiles_);
    for (unsigned int p = 0; p < nbProfiles_; p++) {
      VVdouble * dpxy__son_p = & (*dpxy__son)[p];
      dpxy__son_p->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++) {
	(*dpxy__son_p)[x].resize(nbStates_);
      }
    }

    VVVdouble* d2pxy__son = & d2pxy_[son->getId()];
    d2pxy__son->resize(nbProfiles_);
    for (unsigned int p = 0; p < nbProfiles_; p++) {
      VVdouble * d2pxy__son_p = & (*d2pxy__son)[p];
      d2pxy__son_p->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++) {
	(*d2pxy__son_p)[x].resize(nbStates_);
      }
    }
  }

  if (initialized_) {
    initParameters();
    computeAllTransitionProbabilities();
    fireParameterChanged(getParameters());
  }

}

void RNonHomogeneousMultiProfileTreeLikelihood::initParameters() {

  resetParameters_();

  // Branch lengths
  initBranchLengthsParameters();
  addParameters_(brLenParameters_);

  // Substitution models
  for (auto p : profileSet_) {
    addParameters_(p->getIndependentParameters());
  }
  
  // Profile frequencies
  //for (size_t i=1; i < profileSet_.size(); i++) {
    //addParameter_(profileFreqs_->getParameter("theta" + to_string(i)).clone());
    //}
  addParameters_(profileFreqs_->getIndependentParameters());

}

void RNonHomogeneousMultiProfileTreeLikelihood::initBranchLengthsParameters() {

  brLenParameters_.reset();
  double l1 = 0, l2 = 0;
  for (unsigned int i = 0; i < nbNodes_; i++) {
    double d = minimumBrLen_;
    if (!nodes_[i]->hasDistanceToFather()) {
      nodes_[i]->setDistanceToFather(minimumBrLen_);
    } else {
      d = nodes_[i]->getDistanceToFather();
      if (d < minimumBrLen_) {
	d = minimumBrLen_;
      }
      
      if (d > maximumBrLen_) {
	d = maximumBrLen_;
      }
    }

    brLenParameters_.addParameter(Parameter("BrLen" + TextTools::toString(i),
					    d,
					    brLenConstraint_->clone(),
					    true));
  }

}

void RNonHomogeneousMultiProfileTreeLikelihood::initialize() {

  if (initialized_)
    throw
      Exception("RNonHomogeneousMultiProfileTreeLikelihood::initialize(). Object is already initialized.");
  if (!data_)
    throw
      Exception("RNonHomogeneousMultiProfileTreeLikelihood::initialize(). Data not yet set.");

  initParameters();
  initialized_ = true;
  computeAllTransitionProbabilities();
  fireParameterChanged(getParameters());

}

ParameterList RNonHomogeneousMultiProfileTreeLikelihood::getBranchLengthsParameters() const {

  if (!initialized_) throw Exception("RNonHomogeneousMultiProfileTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");

  return brLenParameters_.getCommonParametersWith(getParameters());

}

ParameterList RNonHomogeneousMultiProfileTreeLikelihood::getSubstitutionModelParameters() const {

  if(!initialized_) throw
     Exception("RNonHomogeneousMultiProfileTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");

  ParameterList pl;
  for (auto p : profileSet_) {
    pl.addParameters(p->getParameters().getCommonParametersWith(getParameters()));
  }
  
  return pl;

}

ParameterList RNonHomogeneousMultiProfileTreeLikelihood::getProfileFrequencyParameters() const {

  if (!initialized_) throw Exception ("RNonHomogeneousMultiProfileTreeLikelihood::getProfileFrequencyParameters(). Object is not initialized.");

  return profileFreqs_->getParameters().getCommonParametersWith(getParameters());

}


// Here we have to sum over the root frequencies of each profile to get the required vector
const vector<double>&  RNonHomogeneousMultiProfileTreeLikelihood::getRootFrequencies(size_t siteIndex) const {

  vector<double> * rf = new vector<double>(nbStates_, 0.0);
  
  for (size_t p=0; p<nbProfiles_; ++p) {
    for (size_t s=0; s<nbStates_; ++s) {
      (*rf)[s] += rootFreqs_[p][s] * profileFreqs_->getProbability(p);
    }
  }

  if (fabs(VectorTools::sum(*rf) - 1.0) >= NumConstants::TINY())
    throw
      Exception("RNonHomogeneousMultiProfileTreeLikelihood::getRootFrequencies.Obtained mixed root frequencies do not sum to one.");
  
  return *rf;

}

void RNonHomogeneousMultiProfileTreeLikelihood::applyParameters() {
  if (!initialized_) throw
		       Exception("RNonHomogeneousMultiProfileTreeLikelihood::applyParameters(). Object not initialized.");

  // Branch lengths
  for (unsigned int i = 0; i < nbNodes_; i++) {
    int id = nodes_[i]->getId();
    const Parameter* brLen = &getParameter(string("BrLen") + TextTools::toString(i));
    if (brLen) nodes_[i]->setDistanceToFather(brLen->getValue());
  }

  // Substitution model
  for (auto p : profileSet_) {
    p->matchParametersValues(getParameters());
  }

  profileFreqs_->matchParametersValues(getParameters());
  
}

void RNonHomogeneousMultiProfileTreeLikelihood::setMinimumBranchLength(double minimum) {
  if (minimum > maximumBrLen_)
    throw Exception("AbstractNonHomogeneousTreeLikelihood::setMinimumBranchLength. Minimum branch length sould be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
  minimumBrLen_ = minimum;
  if (brLenConstraint_.get()) brLenConstraint_.release();
  brLenConstraint_.reset(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true));
  initBranchLengthsParameters();
}

void RNonHomogeneousMultiProfileTreeLikelihood::setMaximumBranchLength(double maximum) {
  if (maximum < minimumBrLen_)
    throw Exception("AbstractNonHomogeneousTreeLikelihood::setMaximumBranchLength. Maximum branch length sould be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
  maximumBrLen_ = maximum;
  if (brLenConstraint_.get()) brLenConstraint_.release();
  brLenConstraint_.reset(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true));
  initBranchLengthsParameters();
}

void RNonHomogeneousMultiProfileTreeLikelihood::computeTreeLikelihood() {
  computeSubtreeLikelihood(tree_->getRootNode());
}

 void RNonHomogeneousMultiProfileTreeLikelihood::computeSubtreeLikelihood(const Node * node) {
   vector <size_t> allProfiles;
   for (size_t i=0; i < nbProfiles_; i++) {
     allProfiles.push_back(i);
   }
   computeSubtreeLikelihoodForProfiles(node, allProfiles); 
}

void RNonHomogeneousMultiProfileTreeLikelihood::computeTreeLikelihoodForProfiles(vector <size_t>& dirtyProfiles) {
  computeSubtreeLikelihoodForProfiles(tree_->getRootNode(), dirtyProfiles);
}

void RNonHomogeneousMultiProfileTreeLikelihood::computeSubtreeLikelihoodForProfiles(const Node* node, vector <size_t>& dirtyProfiles) {

  if (node->isLeaf()) return;

  size_t nbSites  = likelihoodData_->getLikelihoodArray(node->getId()).size();
  size_t nbNodes  = node->getNumberOfSons();

  // Must reset the likelihood array first (i.e. set all of them to 1):
  VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
  for (size_t i = 0; i < nbSites; i++)
  {
    //For each site in the sequence,
    VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
    for (size_t p : dirtyProfiles)
    {
      //For each profile
      Vdouble* _likelihoods_node_i_p = &(*_likelihoods_node_i)[p];
      for (size_t x = 0; x < nbStates_; x++)
      {
        //For each initial state,
        (*_likelihoods_node_i_p)[x] = 1.;
      }
    }
  }

  for (size_t l = 0; l < nbNodes; l++)
  {
    //For each son node,

    const Node* son = node->getSon(l);

    computeSubtreeLikelihoodForProfiles(son, dirtyProfiles); //Recursive method:

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    for (size_t i = 0; i < nbSites; i++)
    {
      //For each site in the sequence,
      VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
      VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
      for (size_t p : dirtyProfiles)
      {
        //For each profile
        Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
        Vdouble* _likelihoods_node_i_p = &(*_likelihoods_node_i)[p];
        VVdouble* pxy__son_p = &(*pxy__son)[p];
        for (size_t x = 0; x < nbStates_; x++)
        {
          //For each initial state,
          Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
          {
            likelihood += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
          }
          (*_likelihoods_node_i_p)[x] *= likelihood;
        }
      }
    }
  }
  
}

//identifies the set of nodes for which parameters need to be recalculated
void RNonHomogeneousMultiProfileTreeLikelihood::fireParameterChanged(const ParameterList& params) {

  applyParameters();

  vector<size_t> dirtyProfiles;

  vector <string> changedBLens = params.getCommonParametersWith(brLenParameters_).getParameterNames();
  vector <int> changedLenNodeIds ;
  for (size_t i = 0; i < changedBLens.size(); i++) {
    int changedId = TextTools::to < size_t> (changedBLens[i].substr(5));
    changedLenNodeIds.push_back(changedId);
    const Node * changedLenNode = nodes_[changedId];
    computeTransitionProbabilitiesForNode(changedLenNode);
  }
  
  for (size_t p=0; p < nbProfiles_ ; p++) {

    vector <int> idsForProfile;
    vector <const Node*> nodesForProfile;
    
    vector <string> tmp =
      params.getCommonParametersWith(profileSet_[p]->getNodeParameters()).getParameterNames();
    for (string i : tmp) {
      vector <int> tmpv = profileSet_[p]->getNodesWithParameter(i);
      idsForProfile = VectorTools::vectorUnion(idsForProfile, tmpv);
    }
    
    for (size_t i = 0; i < idsForProfile.size(); i++) {
      if (find(changedLenNodeIds.begin(), changedLenNodeIds.end(), idsForProfile[i])
	  == changedLenNodeIds.end())
	nodesForProfile.push_back(idToNode_[idsForProfile[i]]);
    }

    for (size_t i = 0; i < nodesForProfile.size(); i++) {
      computeTransitionProbabilitiesForNodeForProfile(nodesForProfile[i], p);
    }

    if (idsForProfile.size() > 0) {
      dirtyProfiles.push_back(p);
    }

    rootFreqs_[p] = profileSet_[p]->getRootFrequencies();

  }
  
  if (changedLenNodeIds.size() > 0) {
    dirtyProfiles.clear();
    for (size_t i=0; i < nbProfiles_; i++) {
      dirtyProfiles.push_back(i);
    }
  }
  
  computeTreeLikelihoodForProfiles(dirtyProfiles);
  
  minusLogLik_ = -getLogLikelihood();

}


void RNonHomogeneousMultiProfileTreeLikelihood::computeTransitionProbabilitiesForNodeForProfile (const Node * node, size_t profile) {

  TransitionModel * modelForNode = profileSet_[profile]->getModelForNode(node->getId()); 

  double l = node->getDistanceToFather();
  RowMatrix <double> Q = modelForNode->getPij_t(l);

  VVdouble * pxy__node_p = &(pxy_[node->getId()][profile]);
  
  for (size_t x=0; x < nbStates_; x++) {
    for (size_t y=0; y < nbStates_; y++) {
      (*pxy__node_p)[x][y] = Q(x, y);
    }
  }

  if (computeFirstOrderDerivatives_) {
    
    VVdouble * dpxy__node_p = &(dpxy_[node->getId()][profile]);
    RowMatrix<double> dQ = modelForNode->getdPij_dt(l);

    for (size_t x=0; x < nbStates_; x++) {
      for (size_t y=0; y < nbStates_; y++) {
	(*dpxy__node_p)[x][y] = dQ(x, y);
      }
    }
    
  }

  if (computeSecondOrderDerivatives_) {

    VVdouble * d2pxy__node_p = &(d2pxy_[node->getId()][profile]);
    RowMatrix<double> d2Q = modelForNode->getd2Pij_dt2(l);

    
    for (size_t x=0; x < nbStates_; x++) {
      for (size_t y=0; y < nbStates_; y++) {
	(*d2pxy__node_p)[x][y] = d2Q(x, y);
      }
    }

  }

}

// Deletes data, imports new data from sites object,  
void RNonHomogeneousMultiProfileTreeLikelihood::setData(const SiteContainer& sites) {

  if (data_) delete data_;

  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());

  likelihoodData_->initLikelihoods(*data_, *(profileSet_[0]->getModelForNode(0)));

  nbSites_ = likelihoodData_->getNumberOfSites();
  nbStates_ = likelihoodData_->getNumberOfStates();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();

  initialized_ = false;
}


void RNonHomogeneousMultiProfileTreeLikelihood::computeAllTransitionProbabilities() {
  for(unsigned int l = 0; l < nbNodes_; l++) {
    Node * node = nodes_[l];
    computeTransitionProbabilitiesForNode(node);
  }
  for (int i=0; i < profileSet_.size(); i++) {
    rootFreqs_[i] = profileSet_[i]->getRootFrequencies();
  }
}

vector <TransitionModel*> RNonHomogeneousMultiProfileTreeLikelihood::getModelsForNode(int nodeId) const {

  vector <TransitionModel*> models;
  for (int i=0; i<nbProfiles_; i++) {
    models.push_back(profileSet_[i]->getModelForNode(nodeId));
  }

  return models;
}
 
void RNonHomogeneousMultiProfileTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node) {

  vector <TransitionModel*> modelsForNode = getModelsForNode(node->getId());
  double l = node->getDistanceToFather();

  VVVdouble * pxy__node = &pxy_[node->getId()];
  
  for (unsigned int p = 0; p < nbProfiles_; p++) {
    TransitionModel* model = modelsForNode[p];
    VVdouble * pxy__node_p = & (*pxy__node)[p];
    RowMatrix<double> Q = model->getPij_t(l);
    for (unsigned int x = 0; x < nbStates_; x++) {
      Vdouble * pxy__node_p_x = & (*pxy__node_p)[x];
      for (unsigned int y = 0; y < nbStates_; y++) {
	(* pxy__node_p_x)[y] = Q(x, y);
      }
    }
  }

  if (computeFirstOrderDerivatives_) {
    VVVdouble * dpxy__node = & dpxy_[node->getId()];

    for(unsigned int p = 0; p < nbProfiles_; p++) {
      VVdouble * dpxy__node_p = & (* dpxy__node)[p];
      TransitionModel* model = modelsForNode[p];

      RowMatrix<double> dQ = model->getdPij_dt(l);

      for(unsigned int x = 0; x < nbStates_; x++) {
	Vdouble * dpxy__node_p_x = & (* dpxy__node_p)[x];
	for (unsigned int y = 0; y < nbStates_; y++)
	  (* dpxy__node_p_x)[y] = dQ(x, y);
      }
    }
  }

  if (computeSecondOrderDerivatives_) {
    VVVdouble * d2pxy__node = & d2pxy_[node->getId()];
    for (unsigned int p = 0; p < nbProfiles_; p++) {
      VVdouble * d2pxy__node_p = & (* d2pxy__node)[p];
      TransitionModel* model = modelsForNode[p];

      RowMatrix<double> d2Q = model->getd2Pij_dt2(l);

      for(unsigned int x = 0; x < nbStates_; x++) {
	Vdouble * d2pxy__node_p_x = & (* d2pxy__node_p)[x];
	for(unsigned int y = 0; y < nbStates_; y++) {
	  (* d2pxy__node_p_x)[y] = d2Q(x, y);
	}
      }
    }
  }
  
}

double RNonHomogeneousMultiProfileTreeLikelihood::getLikelihood() const {

   double l = 1.0;
   for (size_t i = 0; i < nbSites_; i++) {
     l *= getLikelihoodForASite(i);
   }
   return l;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getLogLikelihood() const {
  double ll = 0;
  vector<double> la (nbSites_);
  for (size_t i = 0; i < nbSites_; i++) {
    la[i] = getLogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  for (size_t i = nbSites_; i > 0; i--) {
    ll += la[i - 1];
  }
  return ll;
}
 
double RNonHomogeneousMultiProfileTreeLikelihood::getLikelihoodForASite(size_t site) const {
  double l = 0;
  for (size_t p = 0; p < nbProfiles_; p++) {
    l += getLikelihoodForASiteForAProfile(site, p) * profileFreqs_->getProbability(p);
  }
  return l;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getLogLikelihoodForASite(size_t site) const {
  double l = 0;
  for (size_t p=0; p < nbProfiles_; p++) {
    l += getLikelihoodForASiteForAProfile(site, p) * profileFreqs_->getProbability(p);
  }
  if (l < 0) l = 0;
  return log(l);
}


double RNonHomogeneousMultiProfileTreeLikelihood::getLikelihoodForASiteForAProfileForAState (size_t site, size_t profile, int state) const {

   double l = 0;

   Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][profile];

   l = (*la)[state] * rootFreqs_[profile][state];

   return l;
 }

double RNonHomogeneousMultiProfileTreeLikelihood::getLikelihoodForASiteForAState (size_t site,
										   int state) const {
  //sum across profiles
  double l = 0;
  for (size_t p = 0; p < nbProfiles_; p++) {
    l += getLikelihoodForASiteForAProfileForAState(site, p, state);
  }
  if (l < 0) l = 0;
  return l;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getLikelihoodForASiteForAProfile(size_t site, size_t profile) const {
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][profile];
  for (size_t i=0; i < nbStates_; i++) {
    l += (*la)[i] * rootFreqs_[profile][i];
  }
  return l;
}

 
double RNonHomogeneousMultiProfileTreeLikelihood::getLogLikelihoodForASiteForAProfile(size_t site, size_t profile) const {
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][profile];
  for (size_t i=0; i < nbStates_; i++) {
    l += (*la)[i] * rootFreqs_[profile][i];
  }
  return log(l);
}

 
double RNonHomogeneousMultiProfileTreeLikelihood::getLogLikelihoodForASiteForAState (size_t site, int state) const {

  size_t nbProfiles = nbProfiles_;
  double l = 0;
  for (size_t i=0; i<nbProfiles; ++i) {
    l += getLikelihoodForASiteForAProfileForAState(site, i, state) * profileFreqs_->getProbability(i);
  }
  return log(l);
}

/**
 * !TODO Total hack; no right way to do this. Not for actual use.
 * Will need to rewrite top-level TreeLikelihood classes, joy...
 */
const TransitionModel* RNonHomogeneousMultiProfileTreeLikelihood::getModelForSite(int nodeId, size_t siteIndex) const {

  vector <SubstitutionModel*> modelsForSite;
  for (auto p : profileSet_) {
    modelsForSite.push_back(p->getSubstitutionModelForNode(nodeId));
  }

  Vdouble* vproba = new Vdouble(profileFreqs_->getProbabilities());

  Vdouble* vrate = new Vdouble(4, 1.0);
   
  TransitionModel* branchSpecificMixture = new MixtureOfSubstitutionModels( getAlphabet(), modelsForSite, *vproba, *vrate);

  return branchSpecificMixture;
}
 
/**
 * !TODO: As above.
 */
TransitionModel* RNonHomogeneousMultiProfileTreeLikelihood::getModelForSite(int nodeId, size_t siteIndex) {

  vector <SubstitutionModel*> modelsForSite;
  for (auto p : profileSet_)  {
    modelsForSite.push_back(p->getSubstitutionModelForNode(nodeId));
  }

  Vdouble vproba(profileFreqs_->getProbabilities());

  Vdouble vrate(4, 1.0);
   
  TransitionModel* branchSpecificMixture = new MixtureOfSubstitutionModels( getAlphabet(), modelsForSite, vproba, vrate);

  return branchSpecificMixture;
}


/**
 * !TODO: As above.
 */
VVdouble RNonHomogeneousMultiProfileTreeLikelihood::getTransitionProbabilities(int nodeId, size_t siteIndex) const {

  VVVdouble p3 = getTransitionProbabilitiesPerProfile(nodeId, siteIndex);
  VVdouble p2;
  Vdouble probs = profileFreqs_->getProbabilities();
  p2.resize(getNumberOfStates());
  for (size_t i=0; i<p2.size(); ++i) {
    p2[i].resize(getNumberOfStates());
    for (size_t j=0; j<p2.size(); ++j) {
      for (size_t k=0; k<nbProfiles_; ++k) {
	p2[i][j] += p3[k][i][j] * probs[k];
      }
    }
  }
  return p2;
}
 
ParameterList RNonHomogeneousMultiProfileTreeLikelihood::getDerivableParameters() const {

  if (!initialized_)
    throw
      Exception("AbstractDiscreteRatesAcrossSitesTreeLikelihood::getDerivableParameters(). Object is not initialized.");
    
  return getBranchLengthsParameters();

}

  
ParameterList RNonHomogeneousMultiProfileTreeLikelihood::getNonDerivableParameters() const {

  if (!initialized_)
    throw
      Exception("AbstractDiscreteRatesAcrossSitesTreeLikelihood::getNonDerivableParameters(). Object is not initialized.");
  ParameterList tmp = getSubstitutionModelParameters();
  tmp.addParameters(getProfileFrequencyParameters());

  return tmp;

}


void RNonHomogeneousMultiProfileTreeLikelihood::setParameters(const ParameterList& parameters) {
  setParametersValues(parameters);
}


double RNonHomogeneousMultiProfileTreeLikelihood::getValue() const {

  if (!isInitialized())
    throw
      Exception("RNonHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getFirstOrderDerivative(const string &variable) const {
  if (!hasParameter(variable))
    throw
      ParameterNotFoundException("RNonHomogeneousMultiprofileTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getProfileFrequencyParameters().hasParameter(variable)) {
    throw Exception("Derivatives respective to profile frequency parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable)) {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  const_cast < RNonHomogeneousMultiProfileTreeLikelihood* > (this)->computeTreeDLikelihood(variable);

  return -getDLogLikelihood();
}

double RNonHomogeneousMultiProfileTreeLikelihood::getSecondOrderDerivative(const string &variable) const {
  if (!hasParameter(variable))
    throw
      ParameterNotFoundException("RNonHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getProfileFrequencyParameters().hasParameter(variable)) {
    throw
      Exception("Derivatives respective to profile frequency parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable)) {
    throw Exception("Derivatives respective to substution model parameters are not implemented.");
  }

  (const_cast <RNonHomogeneousMultiProfileTreeLikelihood*> (this))->computeTreeD2Likelihood(variable);

  return -getD2LogLikelihood();
}

 
/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/

double RNonHomogeneousMultiProfileTreeLikelihood::getDLikelihoodForASiteForAProfile(size_t site, size_t profile) const {

  double dl = 0;
  Vdouble* dla = &likelihoodData_->getDLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][profile];
  for (size_t i=0; i<nbStates_; i++) {
    dl += (*dla)[i] * rootFreqs_[profile][i];
  }
  return dl;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getDLikelihoodForASite(size_t site) const {

  double dl=0;
  for (size_t i=0; i < nbProfiles_; i++) {
    dl += getDLikelihoodForASiteForAProfile(site, i) * profileFreqs_->getProbability(i);
  }
  return dl;
}

double RNonHomogeneousMultiProfileTreeLikelihood::getDLogLikelihoodForASite(size_t site) const {

  return getDLikelihoodForASite(site) / getLikelihoodForASite(site); 
  
}

double RNonHomogeneousMultiProfileTreeLikelihood::getDLogLikelihood() const {

  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++) {
    dl += getDLogLikelihoodForASite(i);
  }
  return dl;
}

void RNonHomogeneousMultiProfileTreeLikelihood::computeTreeDLikelihood(const string& variable) {

  // Does different things depending on the type of parameter
  
  if (variable == "RootPosition")
  {
    const Node* father = tree_->getRootNode();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
    size_t nbSites  = _dLikelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++)
    {
      VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
      for (size_t p = 0; p < nbProfiles_; p++)
      {
        Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*_dLikelihoods_father_i_p)[s] = 1.;
        }
      }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++)
    {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double len = getParameterValue("BrLenRoot");

        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoodsroot1__i_p = &(*_likelihoodsroot1__i)[p];
            Vdouble* _likelihoodsroot2__i_p = &(*_likelihoodsroot2__i)[p];
            Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
            VVdouble* dpxy_root1__p  = &(*dpxy_root1_)[p];
            VVdouble* dpxy_root2__p  = &(*dpxy_root2_)[p];
            VVdouble* pxy_root1__p   = &(*pxy_root1_)[p];
            VVdouble* pxy_root2__p   = &(*pxy_root2_)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* dpxy_root1__p_x  = &(*dpxy_root1__p)[x];
              Vdouble* dpxy_root2__p_x  = &(*dpxy_root2__p)[x];
              Vdouble* pxy_root1__p_x   = &(*pxy_root1__p)[x];
              Vdouble* pxy_root2__p_x   = &(*pxy_root2__p)[x];
	      // After all that, actual calculations...
              double dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl1  += (*dpxy_root1__p_x)[y]  * (*_likelihoodsroot1__i_p)[y]; // Summing derivatives across states
                dl2  += (*dpxy_root2__p_x)[y]  * (*_likelihoodsroot2__i_p)[y];
                l1   += (*pxy_root1__p_x)[y]   * (*_likelihoodsroot1__i_p)[y];
                l2   += (*pxy_root2__p_x)[y]   * (*_likelihoodsroot2__i_p)[y];
              }
              double dl = len * (dl1 * l2 - dl2 * l1);
              (*_dLikelihoods_father_i_p)[x] *= dl;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
          VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
            Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
            VVdouble* pxy__son_p = &(*pxy__son)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              double dl = 0;
              Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
              }
              (*_dLikelihoods_father_i_p)[x] *= dl;
            }
          }
        }
      }
    }
    return;
  }

  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t p = 0; p < nbProfiles_; p++)
    {
      Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_p)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* dpxy__son = &dpxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
          VVdouble* dpxy__son_p = &(*dpxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* dpxy__son_p_x = &(*dpxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*dpxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_dLikelihoods_father_i_p)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_dLikelihoods_father_i_p)[x] *= dl;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood(father);

}


void RNonHomogeneousMultiProfileTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t p = 0; p < nbProfiles_; p++)
    {
      Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_p)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _dLikelihoods_son = &likelihoodData_->getDLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _dLikelihoods_son_i = &(*_dLikelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _dLikelihoods_son_i_p = &(*_dLikelihoods_son_i)[p];
          Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_p_x)[y] * (*_dLikelihoods_son_i_p)[y];
            }
            (*_dLikelihoods_father_i_p)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _dLikelihoods_father_i_p = &(*_dLikelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_dLikelihoods_father_i_p)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeDLikelihood(father);
}




/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/

double RNonHomogeneousMultiProfileTreeLikelihood::getD2LikelihoodForASiteForAProfile(
  size_t site,
  size_t profile) const
{
  double d2l = 0;
  Vdouble* d2la = &likelihoodData_->getD2LikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][profile];
  for (size_t i = 0; i < nbStates_; i++)
  {
    d2l += (*d2la)[i] * rootFreqs_[profile][i];
  }
  return d2l;
}

/******************************************************************************/

double RNonHomogeneousMultiProfileTreeLikelihood::getD2LikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t i = 0; i < nbProfiles_; i++)
  {
    d2l += getD2LikelihoodForASiteForAProfile(site, i) * profileFreqs_->getProbability(i);
  }
  return d2l;
}

/******************************************************************************/

double RNonHomogeneousMultiProfileTreeLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/

double RNonHomogeneousMultiProfileTreeLikelihood::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++)
  {
    dl += getD2LogLikelihoodForASite(i);
  }
  return dl;
}

/******************************************************************************/


void RNonHomogeneousMultiProfileTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  if (variable == "BrLenRoot")
  {
    const Node* father = tree_->getRootNode();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
    size_t nbSites  = _d2Likelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++)
    {
      VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
      for (size_t p = 0; p < nbProfiles_; p++)
      {
        Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*_d2Likelihoods_father_i_p)[s] = 1.;
        }
      }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++)
    {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double pos = getParameterValue("RootPosition");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoodsroot1__i_p = &(*_likelihoodsroot1__i)[p];
            Vdouble* _likelihoodsroot2__i_p = &(*_likelihoodsroot2__i)[p];
            Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
            VVdouble* d2pxy_root1__p = &(*d2pxy_root1_)[p];
            VVdouble* d2pxy_root2__p = &(*d2pxy_root2_)[p];
            VVdouble* dpxy_root1__p  = &(*dpxy_root1_)[p];
            VVdouble* dpxy_root2__p  = &(*dpxy_root2_)[p];
            VVdouble* pxy_root1__p   = &(*pxy_root1_)[p];
            VVdouble* pxy_root2__p   = &(*pxy_root2_)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__p_x = &(*d2pxy_root1__p)[x];
              Vdouble* d2pxy_root2__p_x = &(*d2pxy_root2__p)[x];
              Vdouble* dpxy_root1__p_x  = &(*dpxy_root1__p)[x];
              Vdouble* dpxy_root2__p_x  = &(*dpxy_root2__p)[x];
              Vdouble* pxy_root1__p_x   = &(*pxy_root1__p)[x];
              Vdouble* pxy_root2__p_x   = &(*pxy_root2__p)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l1 += (*d2pxy_root1__p_x)[y] * (*_likelihoodsroot1__i_p)[y];
                d2l2 += (*d2pxy_root2__p_x)[y] * (*_likelihoodsroot2__i_p)[y];
                dl1  += (*dpxy_root1__p_x)[y]  * (*_likelihoodsroot1__i_p)[y];
                dl2  += (*dpxy_root2__p_x)[y]  * (*_likelihoodsroot2__i_p)[y];
                l1   += (*pxy_root1__p_x)[y]   * (*_likelihoodsroot1__i_p)[y];
                l2   += (*pxy_root2__p_x)[y]   * (*_likelihoodsroot2__i_p)[y];
              }
              double d2l = pos * pos * d2l1 * l2 + (1. - pos) * (1. - pos) * d2l2 * l1 + 2 * pos * (1. - pos) * dl1 * dl2;
              (*_d2Likelihoods_father_i_p)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
            Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
            VVdouble* pxy__son_p = &(*pxy__son)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              double d2l = 0;
              Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
              }
              (*_d2Likelihoods_father_i_p)[x] *= d2l;
            }
          }
        }
      }
    }
    return;
  }
  else if (variable == "RootPosition")
  {
    const Node* father = tree_->getRootNode();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
    size_t nbSites  = _d2Likelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++)
    {
      VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
      for (size_t p = 0; p < nbProfiles_; p++)
      {
        Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*_d2Likelihoods_father_i_p)[s] = 1.;
        }
      }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++)
    {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double len = getParameterValue("BrLenRoot");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoodsroot1__i_p = &(*_likelihoodsroot1__i)[p];
            Vdouble* _likelihoodsroot2__i_p = &(*_likelihoodsroot2__i)[p];
            Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
            VVdouble* d2pxy_root1__p = &(*d2pxy_root1_)[p];
            VVdouble* d2pxy_root2__p = &(*d2pxy_root2_)[p];
            VVdouble* dpxy_root1__p  = &(*dpxy_root1_)[p];
            VVdouble* dpxy_root2__p  = &(*dpxy_root2_)[p];
            VVdouble* pxy_root1__p   = &(*pxy_root1_)[p];
            VVdouble* pxy_root2__p   = &(*pxy_root2_)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__p_x = &(*d2pxy_root1__p)[x];
              Vdouble* d2pxy_root2__p_x = &(*d2pxy_root2__p)[x];
              Vdouble* dpxy_root1__p_x  = &(*dpxy_root1__p)[x];
              Vdouble* dpxy_root2__p_x  = &(*dpxy_root2__p)[x];
              Vdouble* pxy_root1__p_x   = &(*pxy_root1__p)[x];
              Vdouble* pxy_root2__p_x   = &(*pxy_root2__p)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l1 += (*d2pxy_root1__p_x)[y] * (*_likelihoodsroot1__i_p)[y];
                d2l2 += (*d2pxy_root2__p_x)[y] * (*_likelihoodsroot2__i_p)[y];
                dl1  += (*dpxy_root1__p_x)[y]  * (*_likelihoodsroot1__i_p)[y];
                dl2  += (*dpxy_root2__p_x)[y]  * (*_likelihoodsroot2__i_p)[y];
                l1   += (*pxy_root1__p_x)[y]   * (*_likelihoodsroot1__i_p)[y];
                l2   += (*pxy_root2__p_x)[y]   * (*_likelihoodsroot2__i_p)[y];
              }
              double d2l = len * len * (d2l1 * l2 + d2l2 * l1 - 2 * dl1 * dl2);
              (*_d2Likelihoods_father_i_p)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t p = 0; p < nbProfiles_; p++)
          {
            Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
            Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
            VVdouble* pxy__son_p = &(*pxy__son)[p];
            for (size_t x = 0; x < nbStates_; x++)
            {
              double d2l = 0;
              Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
              }
              (*_d2Likelihoods_father_i_p)[x] *= d2l;
            }
          }
        }
      }
    }
    return;
  }

  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t p = 0; p < nbProfiles_; p++)
    {
      Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_p)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* d2pxy__son = &d2pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
          VVdouble* d2pxy__son_p = &(*d2pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* d2pxy__son_p_x = &(*d2pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*d2pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_d2Likelihoods_father_i_p)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_d2Likelihoods_father_i_p)[x] *= d2l;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void RNonHomogeneousMultiProfileTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t p = 0; p < nbProfiles_; p++)
    {
      Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_p)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _d2Likelihoods_son = &likelihoodData_->getD2LikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _d2Likelihoods_son_i = &(*_d2Likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _d2Likelihoods_son_i_p = &(*_d2Likelihoods_son_i)[p];
          Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_p_x)[y] * (*_d2Likelihoods_son_i_p)[y];
            }
            (*_d2Likelihoods_father_i_p)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t p = 0; p < nbProfiles_; p++)
        {
          Vdouble* _likelihoods_son_i_p = &(*_likelihoods_son_i)[p];
          Vdouble* _d2Likelihoods_father_i_p = &(*_d2Likelihoods_father_i)[p];
          VVdouble* pxy__son_p = &(*pxy__son)[p];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_p_x = &(*pxy__son_p)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_p_x)[y] * (*_likelihoods_son_i_p)[y];
            }
            (*_d2Likelihoods_father_i_p)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/
