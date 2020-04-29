
/** File: BranchMixtureMutselInference.cpp **/

/** Contains all inference, diagnostic and I/O procedures for branch-specific mixture non-equilibrium  ** mutation-selection modelling. **/

#include <numeric>
#include <boost/math/distributions/chi_squared.hpp>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/CodonSiteTools.h>

#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/MixtureOfSubstitutionModels.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

#include "MutationSelectionModel.h"
#include "AdditionalOptimizationTools.h"
#include "BranchMixtureMutSelInference.h"

using namespace std;
using namespace bpp;

/**********************************************************************************************/

// Constructors, destructors, operators

BranchMixtureMutSelInference::BranchMixtureMutSelInference(const Tree * tree,
							   const SiteContainer * sites,
							   const vector<int>& assumedBreakpoints):
  tree_(tree),
  sites_(sites),
  assumedBreakpoints_(assumedBreakpoints),
  inferredBreakpoints_()
{

  init_();

}

BranchMixtureMutSelInference::BranchMixtureMutSelInference(const BranchMixtureMutSelInference& inf):
  tree_(inf.tree_),
  sites_(inf.sites_),
  maxBreakpoints_(inf.maxBreakpoints_),
  assumedBreakpoints_(inf.assumedBreakpoints_),
  inferredBreakpoints_(inf.inferredBreakpoints_),
  nhtl_(inf.nhtl_),
  branchIds_(inf.branchIds_),
  positionLikelihoods_(inf.positionLikelihoods_),
  aiccs_(inf.aiccs_)
{}

BranchMixtureMutSelInference& BranchMixtureMutSelInference::operator= (const BranchMixtureMutSelInference& inf) {
  tree_ = inf.tree_;
  sites_ = inf.sites_;
  maxBreakpoints_ = inf.maxBreakpoints_;
  assumedBreakpoints_ = inf.assumedBreakpoints_;
  inferredBreakpoints_ = inf.inferredBreakpoints_;
  nhtl_ = inf.nhtl_->clone();
  branchIds_ = inf.branchIds_;
  positionLikelihoods_ = inf.positionLikelihoods_;
  aiccs_ = inf.aiccs_;

  return *this;

}

BranchMixtureMutSelInference::~BranchMixtureMutSelInference() {

  delete nhtl_;
  delete storedLikelihood_;

}

/**********************************************************************************************/

// Main inference procedures


/** Deconvolute a mixture of mutsel models over sites with an unknown number of components.
 ** Largely intractable for any more than one component (two max). With 1 initial model and
 ** nsRootFreqs=true, this estimates codon frequencies at the root for a non-stationary model.
 **/
void BranchMixtureMutSelInference::bppInferBackgroundMixture(size_t maxInitialModels,
							     bool nsRootFreqs) {

  if (nsRootFreqs) {

    CodonAlphabet * codonAlpha = new CodonAlphabet(&AlphabetTools::DNA_ALPHABET);
    StandardGeneticCode * gCode = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
    FullCodonFrequenciesSet * freqs = new FullCodonFrequenciesSet(gCode);
    
    SubstitutionModelSet * newSet = new SubstitutionModelSet(codonAlpha, freqs);
    SubstitutionModelSet * oldSet = nhtl_->getSubstitutionModelSet();

    for (int i=0; i < oldSet->getNumberOfModels(); i++) {

      newSet->addModel(oldSet->getModel(i)->clone(), oldSet->getNodesWithModel(i));
      
    }
    
    nhtl_->setSubstitutionModelSet(newSet);

  }
  
  optimizeAllParameters_();
  
  optimizeBPMixture_(tree_->getRootId(), maxInitialModels, true);

}

/** Assume a mixture with identical components over time, but where some sites have moved between
 ** components (selective regimes). Iteratively find those positions by AICc. The (known) parameters 
 ** of the
 ** mixture components must currently be passed in as a Bio++ parameter list. Experimental.
 **/

void BranchMixtureMutSelInference::bppFindSiteFrequencyShifts(const ParameterList& initialParams,
							      size_t nBackgroundCats,
							      size_t maxFreqShifts) {


  addModelsToBp_(tree_->getRootId(), nBackgroundCats - 1);
  
  nhtl_->matchParametersValues(initialParams);

  double startLnL = optimizeSiteFrequenciesAllBreakpoints_();
  
  double aicc0 = aicc_corrected_(startLnL, getCurrentDimension_(), sampleSize_);
  double aiccPrev = aicc0;
  double aiccCurrent = aicc0;
  aiccs_.insert(pair<int, double>(-1, aicc0));

  double currentBreakpoints = assumedBreakpoints_.size();

  if (currentBreakpoints >= maxFreqShifts)
    throw
      Exception("BranchMixtureMutSelInference::bppFindSiteFrequencyShifts(). Number of initially assumed shifts is  greater than or equal to the  provided maximum breakpoint cap.");

  while (currentBreakpoints < maxFreqShifts) {

    storeLikelihood_();
    
    int propShift = optimizeBreakpointPosition_(0, false);
    double lnLCurrent = optimizeSiteFrequenciesAllBreakpoints_();
    aiccCurrent = aicc_corrected_(lnLCurrent, getCurrentDimension_(), sampleSize_);
    
    aiccs_.insert(pair<int, double>(propShift, aiccCurrent));
    cout << "Trying breakpoint at position " << propShift << ": AICc = " << aiccCurrent << " - ";
    
    if (aiccCurrent < aiccPrev) {
      inferredBreakpoints_.push_back(propShift);
      aiccPrev = aiccCurrent;
      currentBreakpoints++;
      cout << "accepted." << endl;
    } else {
      revertLikelihood_();
      cout << "rejected." << endl;
      break;
    }
    
  }
 
}

/** Test for the presence of a breakpoint at a specific branch, with full optimisation of
 ** parameters. Untested.
 **/

double BranchMixtureMutSelInference::bppTestSpecificBranch(int testBranchIndex,
							 size_t maxInitialModels,
							 size_t maxBreakpointModels) {

  /** optimize initial mixture components
   ** find aicc
   ** for the test branch:
   **     add a new model
   **     optimize
   **     record likelihood
   **     repeat
   ** optimize test branch model
   ** LRT vs initial mixture
   **/

  double startLnl = optimizeBPMixture_(tree_->getRootId(), maxInitialModels, true);
  int startDim = getCurrentDimension_();
  
  double finalLnl = optimizeBPMixture_(testBranchIndex, maxBreakpointModels, true);
  int endDim = getCurrentDimension_();
  
  return likelihoodRatioTest(finalLnl, startLnl, endDim - startDim); 

}

/** FindSiteFrequencyShifts plus full parameter optimisation; needless to say, intractable.
 ** Call with maxInitialModels = maxBreakpointModels = 1 for old-fashioned estimation of 
 ** breakpoint positions. 
 **/
double BranchMixtureMutSelInference::bppTestBestBranch(size_t maxBreakpoints,
						       size_t maxInitialModels,
						       size_t maxBreakpointModels) {

  /** optimize initial mixture components
   ** find aicc
   ** for each clade:
   **     add a new model
   **     optimize
   **     record likelihood, update max
   **     (repeat?)
   ** for the max:
   **     optimize new model and frequencies
   **     LRT vs  initial mixture
   ** potentially rinse and repeat?
   **/
  
  double startLnL = optimizeBPMixture_(tree_->getRootId(), maxInitialModels, true);
  int startDim = getCurrentDimension_();
  double aicc0 =  aicc_corrected_(startLnL, startDim, sampleSize_);
  double aiccPrev = aicc0;
  double aiccCurrent = aicc0;
  aiccs_.insert(pair<int, double>(-1, aicc0));

  double currentBreakpoints = assumedBreakpoints_.size();
  
  if (currentBreakpoints >= maxBreakpoints)
    throw
      Exception("BranchMixtureMutSelInference::bppTestBestBranch(). Number of initially assumed breakpoints is  greater than or equal to the  provided maximum breakpoint cap.");
  
  do {

    storeLikelihood_();
    
    int propBp = optimizeBreakpointPosition_(0, true);
    double branchLnL = optimizeBPMixture_(propBp, maxBreakpointModels, true);
    double lnLCurrent = nhtl_->getLogLikelihood();
    
    aiccCurrent = aicc_corrected_(lnLCurrent, getCurrentDimension_(), sampleSize_);
    aiccs_.insert(pair<int, double>(propBp, aiccCurrent));
    inferredBreakpoints_.push_back(propBp);
    
    cout << "Trying breakpoint at position " << propBp << ": AICc = " << aiccCurrent << " - ";
    
    if (aiccCurrent < aiccPrev) {
      aiccPrev = aiccCurrent;
      currentBreakpoints++;
      cout << "accepted." << endl;
      
    } else {
      revertLikelihood_();
      cout << "rejected." << endl;
      break;
    }

  } while (currentBreakpoints < maxBreakpoints);

  double finalLnL = nhtl_->getLogLikelihood();
  int endDim = getCurrentDimension_();

  return 1.0;
  
}




/**********************************************************************************************/

// I/O procedures

void BranchMixtureMutSelInference::writeParameterEstimates(string path) {
 
  ParameterList pl2 = nhtl_->getParameters();
  ofstream paramOut;
  
  paramOut.open(path.c_str());
  pl2.printParameters(paramOut);
  paramOut.close();
  
  cout << "Parameter estimates written to " << path << "." << endl;

}


void BranchMixtureMutSelInference::writeBreakpointEstimates(string path) {

  ofstream bpOut;

  bpOut.open(path.c_str());
  bpOut << "Breakpoints assumed at nodes: ";
  for (int i : assumedBreakpoints_) bpOut << " " << i;
  bpOut << endl;

  bpOut << "AICc with only assumed breakpoints: ";
  bpOut << aiccs_[-1];

  bpOut << "Breakpoints tested at nodes: ";
  for (int i=0; i < inferredBreakpoints_.size(); i++) {
    bpOut << " " << inferredBreakpoints_[i];
    bpOut << " AICc " << aiccs_[inferredBreakpoints_[i]];
    bpOut << " - " << (i == inferredBreakpoints_.size() - 1)? "rejected;" : "accepted;";
    bpOut << " LogLikelihoods ";
    for (int l=0; l < positionLikelihoods_[inferredBreakpoints_[i]].size(); l++)
      bpOut << " " <<  TextTools::toString(positionLikelihoods_[inferredBreakpoints_[i]][l]);
    bpOut << endl;
  }
  bpOut << endl;

  bpOut.close();

  cout << "Breakpoints written to " << path << "." << endl; 
  
}


void BranchMixtureMutSelInference::writeEquilibriumFrequencies(string path) {

  ofstream freqOut;
  freqOut.open(path.c_str());
  freqOut << "Model Equilibrium Frequencies" << endl;
  
  const VVdouble rootFreqs;

  vector<SubstitutionModel*> mixcomps;
  getAllMixtureComponents(mixcomps);
  
  for (auto m : mixcomps) {

    string modName = m->getName();
    Vdouble rootFreqs = m->getFrequencies();

    freqOut << modName;
    for (double d : rootFreqs) {
      freqOut << "\t" << d;
    }

    freqOut << endl;
  }

  freqOut.close();

}


void BranchMixtureMutSelInference::writeTimings(string path) {

  ofstream timeOut;

  timeOut.open(path.c_str());
  timeOut << "Step,Time,Duration" << endl;
  time_t last=-1;
  for (auto p : timings_) {
    string dur = last >=0? TextTools::toString(difftime(p.second, last)) : "--";
    timeOut << p.first << "," << p.second << "," << dur << endl;
    last = p.second;
  }

  timeOut.close();

}


void BranchMixtureMutSelInference::setAssumedBreakpoints(const vector<int>& breakpoints) {

  assumedBreakpoints_ = breakpoints;
  resetLikelihood_();

}


void BranchMixtureMutSelInference::getAllMixtureComponents(vector <SubstitutionModel*>& models) const {

  vector <SubstitutionModel*> mixtures;
  vector <string> distinctModels;
  
  for (int i=0; i < branchIds_.size(); i++)
    mixtures.push_back(nhtl_->getSubstitutionModelSet()->getSubstitutionModelForNode(i));

  for (SubstitutionModel * m : mixtures) {

    MixedSubstitutionModel* mix = dynamic_cast<MixedSubstitutionModel*>(m);
    
    for (int i=0; i < mix->getNumberOfModels(); i++) {

      SubstitutionModel * subModel = mix->getNModel(i);

      if ( find(distinctModels.begin(), distinctModels.end(), subModel->getName())
	   == distinctModels.end()) {

	models.push_back(m);
	
      } else {
	distinctModels.push_back(subModel->getName());
      }
    }
    
  }

}


double BranchMixtureMutSelInference::computeLogLikelihood() {

  return nhtl_->getLogLikelihood();

}


double BranchMixtureMutSelInference::computeLogLikelihoodForParameters(const ParameterList& pl,
								       int nbCategories) {

  storeLikelihood_();
  
  vector<int> allBps;
  allBps.push_back(tree_->getRootId());
  for (int bp : assumedBreakpoints_)
    allBps.push_back(bp); 
  for (int bp : inferredBreakpoints_)
    allBps.push_back(bp);
  
  for (int bp : allBps) {

    addModelsToBp_(bp, nbCategories - 1);
    
  }
  
  nhtl_->matchParametersValues(pl);
  
  double paramLike = nhtl_->getLogLikelihood();
  
  revertLikelihood_();

  return paramLike;

}


Vdouble BranchMixtureMutSelInference::computeSiteLogLikelihoodsForParameters(const ParameterList & pl,
									     int nbCategories) {

  storeLikelihood_();
  
  vector<int> allBps;
  allBps.push_back(tree_->getRootId());
  for (int bp : assumedBreakpoints_)
    allBps.push_back(bp); 
  for (int bp : inferredBreakpoints_)
    allBps.push_back(bp);

  for (int bp : allBps) {
  
    addModelsToBp_(bp, nbCategories - 1);
    
  }
  
  nhtl_->matchParametersValues(pl);

  Vdouble paramSiteLikes = nhtl_->getLogLikelihoodForEachSite();

  revertLikelihood_();
  
  return paramSiteLikes;
  
}


void BranchMixtureMutSelInference::getEmpiricalCodonFreqs(VVdouble& tipFreqs,
						       Vdouble& averageFreqs) {

  size_t nSeqs = sites_->getNumberOfSequences();

  for (size_t i=0; i < nhtl_->getNumberOfStates(); i++)
    averageFreqs.push_back(0);
  
  for (size_t i=0; i < nSeqs; i++) {

    map<int, double> codonStateMap;
    
    CodonSiteTools::getFrequencies(sites_->getSequence(i), codonStateMap);

    Vdouble thisTipFreqs;
    for (size_t p=0; p < codonStateMap.size(); p++) {
      thisTipFreqs.push_back(codonStateMap[p]);
      averageFreqs[p] += codonStateMap[p];
    }
    tipFreqs.push_back(thisTipFreqs);
    
  }

  for (size_t i=0; i < nhtl_->getNumberOfStates(); i++) {
    averageFreqs[i] /= nSeqs;
  }

}

/** !TODO: Does not work right as of 03/29/19; must have bumped something. 
 ** 
 ** This method is for diagnostics; it calculates the marginal distribution of codon states
 ** at the tips of the tree under a given model (e.g., true or estimated). 
 **/
void BranchMixtureMutSelInference::computeTipDistrsForParameters(const ParameterList& pl,
								 int nbCategories,
								 VVdouble& outDists) {

  storeLikelihood_();

  vector<int> allBps;
  allBps.push_back(tree_->getRootId());
  for (int bp : assumedBreakpoints_)
    allBps.push_back(bp); 
  for (int bp : inferredBreakpoints_)
    allBps.push_back(bp);
  
  for (int bp : allBps) {

    addModelsToBp_(bp, nbCategories - 1);
    
  }
  
  nhtl_->matchParametersValues(pl);

  Vdouble rootFrequencies = nhtl_->getSubstitutionModelSet()->getRootFrequencies();

  double sum = std::accumulate(rootFrequencies.begin(), rootFrequencies.end(), 0);

  for (int i : tree_->getSonsId(tree_->getRootId())) {
    computeTipDistrs_(i, new Vdouble(rootFrequencies), outDists);
  }

  revertLikelihood_();
  
}


void BranchMixtureMutSelInference::computeTipDistrs_(int branch, Vdouble * startDist, VVdouble& outDists) {

  double branchLength = tree_->getDistanceToFather(branch);

  SubstitutionModel * m = nhtl_->getSubstitutionModelSet()->getSubstitutionModelForNode(branch);
  Vdouble * endDist = new Vdouble(m->getNumberOfStates());

  for (size_t x=0; x < m->getNumberOfStates(); x++) {
    for (size_t y = 0; y < m->getNumberOfStates(); y++)
      (*endDist)[x] += startDist->at(y) * m->Pij_t(y, x, branchLength);
  }

  double sum = std::accumulate(endDist->begin(), endDist->end(), 0.0);

  
  if (tree_->isLeaf(branch)) {
    outDists.push_back(*endDist);
  } else {
    for (int i : tree_->getSonsId(branch)) {
      computeTipDistrs_(i, endDist, outDists);
    }
  }

}


/**********************************************************************************************/

// Private utility procedures

void BranchMixtureMutSelInference::init_() {

  vector<int> ids = tree_->getNodesId();
  int rootId = tree_->getRootId();
  size_t pos = 0;
  for (size_t i=0; i < ids.size(); i++) {
    if (ids[i] == rootId) {
      pos = i;
      break;
    }
  }
  ids.erase(ids.begin() + pos);

  branchIds_ = ids;

  CodonAlphabet * codonAlpha = new CodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  StandardGeneticCode * standardCode = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
  SubstitutionModelSet * modelSet = new SubstitutionModelSet(codonAlpha);

  NucleotideSubstitutionModel* nucModel = new HKY85(&AlphabetTools::DNA_ALPHABET);
  string mixName = "bp" + TextTools::toString(branchIds_.size()) + "_";
  string modName = "m0_";
  SubstitutionModel* addedModel = new MutationSelectionModel(standardCode, nucModel, modName);

  vector<SubstitutionModel*> startModel;
  startModel.push_back(addedModel);
  MixtureOfSubstitutionModels * startMixture = new MixtureOfSubstitutionModels(codonAlpha, startModel);
  startMixture->setNamespace(mixName);
  for (int i=0; i < startMixture->getNumberOfModels(); i++) {
    startMixture->getNModel(0)->setNamespace(mixName +  TextTools::toString(i+1) + "_" + modName);
  }
  modelSet->addModel(startMixture, branchIds_);

  DiscreteDistribution * noRateHet = new ConstantDistribution(1.0);
  
  nhtl_ = new DRNonHomogeneousTreeLikelihood(*tree_,
					     *sites_,
					     modelSet,
					     noRateHet,
					     true,
					     false
					     );
  nhtl_->initialize();

  storedLikelihood_ = nhtl_->clone();
  
}


void BranchMixtureMutSelInference::storeLikelihood_() {
  
  *storedLikelihood_ = *nhtl_;

}


void BranchMixtureMutSelInference::revertLikelihood_() {

  *nhtl_ = *storedLikelihood_;

}


void BranchMixtureMutSelInference::resetLikelihood_() {

  CodonAlphabet * codonAlpha = new CodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  StandardGeneticCode * standardCode = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
  SubstitutionModelSet * modelSet = new SubstitutionModelSet(codonAlpha);
  
  for (int bp=-1; bp < (int)assumedBreakpoints_.size(); bp++) {
      NucleotideSubstitutionModel* nucModel = new HKY85(&AlphabetTools::DNA_ALPHABET);
      int bpnum = bp < 0 ? branchIds_.size() : assumedBreakpoints_[bp]; 
      string mixName = "bp" + TextTools::toString(bpnum) + "_";
      string modName = "m0_";
      SubstitutionModel * addedModel = new MutationSelectionModel(standardCode, nucModel, modName);
  
      vector<SubstitutionModel*> startModel;
      startModel.push_back(addedModel);
      MixtureOfSubstitutionModels * startMixture = new MixtureOfSubstitutionModels(codonAlpha, startModel);
      startMixture->setNamespace(mixName);
      for (int i=0; i < startMixture->getNumberOfModels(); i++) {
	startMixture->getNModel(i)->setNamespace(mixName + TextTools::toString(i+1) + "_" + modName);
      }

      vector <int> addedBranches =
	bp < 0 ? branchIds_ : getDescendantsWithModel(*modelSet, assumedBreakpoints_[bp]);
      
      modelSet->addModel(startMixture, addedBranches);
    
  }


  DiscreteDistribution * noRateHet = new ConstantDistribution(1.0);

  DRNonHomogeneousTreeLikelihood newLike(*tree_,
					 *sites_,
					 modelSet,
					 noRateHet,
					 true,
					 false
					 );
  newLike.initialize();
  
  *nhtl_ = newLike;

  *storedLikelihood_ = *nhtl_; 

}


int BranchMixtureMutSelInference::optimizeBreakpointPosition_(size_t nbAddedModels,
							      bool optimizeParameters) {

  /**
   ** for each clade:
   **     add a new model
   **     optimize
   **     record likelihood, update max
   **     (repeat?)
   **/

  double startLnL = nhtl_->getLogLikelihood();
  int currentBranch = branchIds_.back();
  int currentMLBranch = -1;
  double currentMaxLnL;

  vector <double> lnLs(branchIds_.size());
  VectorTools::fill(lnLs, -NumConstants::VERY_BIG());

  storeLikelihood_();
  
  for (auto it = branchIds_.rbegin(); it != branchIds_.rend(); ++it) {
    
    if (!branchHasBreakpoint_(*it)) {
      
      size_t nbCurrModels =
	dynamic_cast<MixedSubstitutionModel * >(nhtl_->getModelForNode(*it))->getNumberOfModels();

      newBreakpointAt_(*it);
      
      double currentLnL
	= optimizeBPMixture_(*it, nbCurrModels + nbAddedModels, optimizeParameters);
      
      lnLs[*it] = currentLnL;
      if ((currentMLBranch < 0 || currentLnL > currentMaxLnL)) {

	currentMLBranch = *it;
	currentMaxLnL = currentLnL;

      }

      currentBranch = *it;
      
      revertLikelihood_();
      
    }
    
  }

  if (currentMLBranch >= 0)
    newBreakpointAt_(currentMLBranch);

  addModelsToBp_(currentMLBranch, nbAddedModels);

  positionLikelihoods_.insert(pair<int, vector<double> >(currentMLBranch, lnLs));
  
  return currentMLBranch;
  
}


double BranchMixtureMutSelInference::optimizeBPMixture_(int bpPosition,
							size_t maxModels,
							bool optimizeParameters) {

  double startLnL = nhtl_->getLogLikelihood();
  
  double aicc0 = aicc_corrected_(startLnL, getCurrentDimension_(), sampleSize_);
  double aiccPrev = aicc0;
  double aiccCurrent = aicc0;

  cout << "Starting AICc for breakpoint "
       << TextTools::toString(bpPosition)+ ": "
       << TextTools::toString(aicc0)
       << endl << endl;

  int currentModels =
    (bpPosition == tree_->getRootId())?
    dynamic_cast<MixedSubstitutionModel * >(nhtl_->getSubstitutionModelSet()->getModel(0))->getNumberOfModels():
    dynamic_cast<MixedSubstitutionModel * >(nhtl_->getModelForNode(bpPosition))->getNumberOfModels();
  
  while (currentModels < maxModels) {
    
    storeLikelihood_();
    
    addModelsToBp_(bpPosition, 1);
    vector<size_t> optModels;

    if (optimizeParameters) {
      for (int i = 0; i <= currentModels; i++) {
	optModels.push_back(i);
      }
    }
    
    double currentLnL = optimizeMixtureParams_(bpPosition,
					       0.0001,
					       optModels,
					       true);
    aiccCurrent = aicc_corrected_(currentLnL, getCurrentDimension_(), sampleSize_);

    cout << "Trying new background component: AICc = " << aiccCurrent << " - ";
    
    if (aiccCurrent < aiccPrev) {
      aiccPrev = aiccCurrent;
      currentModels++;
      cout << "accepted." << endl << endl;
      
    } else {
      revertLikelihood_();
      cout << "rejected." << endl << endl;
      break;
    }
    
  } 

  vector<size_t> allSubmods;
  for (size_t i=0; i < currentModels; i++)
    allSubmods.push_back(i); 
  double endLnL = optimizeMixtureParams_(bpPosition, 0.00001, allSubmods, true);

  double aiccFinal = aicc_corrected_(endLnL, getCurrentDimension_(), sampleSize_);
  cout << "Final AICc for breakpoint "
       << TextTools::toString(bpPosition) + ": "
       << TextTools::toString(aiccFinal)
       << endl << endl;
  
  return endLnL;

}

/** This is the method that causes all the trouble. It takes a very long time to fully optimise
 ** even a small mixture of mut-sel models.
 **/
double BranchMixtureMutSelInference::optimizeMixtureParams_(int bpPosition,
							    double tolerance,
							    vector <size_t> submodels,
							    bool includeBranchLengths) {

  ParameterList * freqParams = getFrequencyParamsForBP_(bpPosition);
  ParameterList * submodelParams = getParametersForBPSubmodels_(bpPosition, submodels);
  MixedSubstitutionModel * mix =
    (bpPosition == tree_->getRootId())?
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getSubstitutionModelSet()->getModel(0)):
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getModelForNode(bpPosition));

  if (includeBranchLengths) submodelParams->addParameters(nhtl_->getBranchLengthsParameters());
 
  timings_.push_back(pair<string, time_t>("OptMix_" + mix->getName() + "_Init",
					 time(NULL))); 

  submodelParams->addParameters(*freqParams);
  
  optimizeParams_(*submodelParams, 0.000001);

  timings_.push_back(pair<string, time_t>("OptMix_" + mix->getName() + "_Finish",
					 time(NULL)));

 double endLnL = nhtl_->getLogLikelihood();
  
 return endLnL;

}


double BranchMixtureMutSelInference::optimizeAllParameters_() {

  ParameterList allParameters;

  for (int i=0; i < nhtl_->getParameters().size(); i++) {
    if (nhtl_->getParameters()[i].getName().find("relrate") == string::npos)
      allParameters.addParameter(nhtl_->getParameters()[i]);
  }
  
  double endLnL = optimizeParams_(allParameters, 0.00001);
  return endLnL;

}


double BranchMixtureMutSelInference::optimizeParams_(ParameterList optParams, double tolerance) {
 
  timings_.push_back(pair<string, time_t>("OptParams_Init",
					 time(NULL))); 
  
  AdditionalOptimizationTools::optimizeNumericalParametersNHMixed(nhtl_,
								  optParams,
								  0, 1,
								  tolerance,
								  1000000,
								  ApplicationTools::message.get(),
								  ApplicationTools::message.get(),
								  true,
								  1,
								  OptimizationTools::OPTIMIZATION_NEWTON,
								  OptimizationTools::OPTIMIZATION_BRENT);

 timings_.push_back(pair<string, time_t>("OptParams_Finish",
					 time(NULL)));

 double endLnL = nhtl_->getLogLikelihood();
  
 return endLnL;
 
}


double BranchMixtureMutSelInference::optimizeSiteFrequenciesAllBreakpoints_() {

  vector<int> allBreakpoints = VectorTools::vectorUnion(assumedBreakpoints_,
							inferredBreakpoints_);

  ParameterList allFreqParams;

  
  for (int i : allBreakpoints) {
    vector<size_t> dummy;
    allFreqParams.addParameters(*getParametersForBPSubmodels_(i, dummy));
    
  }

  allFreqParams.addParameters(nhtl_->getBranchLengthsParameters());
  
  double endLnL = optimizeParams_(allFreqParams, 0.000001);

  return endLnL;
  
}

ParameterList  * BranchMixtureMutSelInference::getFrequencyParamsForBP_(int bpPosition) {

  MixedSubstitutionModel * optMix = (bpPosition == tree_->getRootId())?
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getSubstitutionModelSet()->getModel(0)):
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getModelForNode(bpPosition));
  
  ParameterList mixParams = optMix->getParameters();
  ParameterList likeParams = nhtl_->getSubstitutionModelParameters();
  ParameterList * outputParams = new ParameterList();
  ParameterList freqParams;
  
  for (int i=0; i < mixParams.size(); i++) {
    if (mixParams[i].getName().find("relproba") != string::npos)
      freqParams.addParameter(mixParams[i]);
  }

  for (int p=0; p < likeParams.size(); p++) {

    string modName = likeParams[p].getName();
    size_t found = modName.rfind("_");
    string subName = modName.substr(0, found);
    
    if (freqParams.hasParameter(subName)) {
      outputParams->addParameter(nhtl_->getParameter(modName));
    }

  }

  return outputParams;
  
}


ParameterList * BranchMixtureMutSelInference::getParametersForBPSubmodels_(int bpPosition,
									   vector<size_t> submodels)  {

  MixedSubstitutionModel * optMix = (bpPosition == tree_->getRootId())?
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getSubstitutionModelSet()->getModel(0)):
    dynamic_cast<MixedSubstitutionModel*>(nhtl_->getModelForNode(bpPosition));
  
  ParameterList likeParams = nhtl_->getSubstitutionModelParameters();
  ParameterList * outputParams = new ParameterList();

  // Parameters for specified submodels
  
  for (int p=0; p < likeParams.size(); p++) {
    
    string parName = likeParams[p].getName();

    size_t isFit = parName.find("fit");
    size_t isHKY = parName.find("HKY85");

    if (isFit != string::npos || isHKY != string::npos) {
    
      size_t firstbreak = parName.find("_", 0);
      string mixName = parName.substr(0, firstbreak) + "_";
      size_t secondbreak = parName.find("_", firstbreak+1);
      string modNumStr = parName.substr(firstbreak+1, secondbreak-(firstbreak+1));
      int modNum = TextTools::toInt(modNumStr);
      auto isTargetModel = find(submodels.begin(), submodels.end(), modNum - 1);
    
      if (mixName == optMix->getNamespace() && isTargetModel != submodels.end()) {
	outputParams->addParameter(nhtl_->getParameter(parName));
      }

    }

  }
  
  return outputParams;

}

/** This should not be this hard **/
void BranchMixtureMutSelInference::addModelsToBp_(int bpPosition, size_t nbAdditionalModels) {

  // Adding to the root model is an annoying special case!

  bool isRoot = tree_->isRoot(bpPosition);
  vector<int> basalBranches = tree_->getSonsId(tree_->getRootId());
  
  if (isRoot &&
      nhtl_->getSubstitutionModelSet()->getModelIndexForNode(basalBranches[0]) !=
      nhtl_->getSubstitutionModelSet()->getModelIndexForNode(basalBranches[1])) {

    throw Exception ("BranchMixtureMutSelInference::addModelsToBp_.Forbidden attempt to overwrite root model.");

  }

  //if (!isRoot && !branchHasBreakpoint_(bpPosition))
  //throw
  //  Exception ("BranchMixtureMutSelInference::addModelsToBp_.No breakpoint at specified branch.");

  string mixName = "bp" + TextTools::toString(bpPosition) + "_";
  StandardGeneticCode * standardCode = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
  CodonAlphabet * codonAlpha = new CodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  SubstitutionModelSet * newSet = new SubstitutionModelSet(codonAlpha);
  for (int i=0; i < nhtl_->getSubstitutionModelSet()->getNumberOfModels(); i++) {
    SubstitutionModel * newMixModel = nhtl_->getSubstitutionModelSet()->getSubstitutionModel(i)->clone();
    newSet->addModel(newMixModel, nhtl_->getSubstitutionModelSet()->getNodesWithModel(i));
  }
  
  MixedSubstitutionModel * oldMixture =
    isRoot?
    dynamic_cast<MixedSubstitutionModel *>(newSet->getModel(0)):
    dynamic_cast<MixedSubstitutionModel *>(newSet->getModelForNode(bpPosition));
  
  vector <SubstitutionModel *> newModels;
  vector <string> allModNames;
  vector <double> relprobas;
  for (int i=0; i < oldMixture->getNumberOfModels(); i++) {
    string modName = "m" + TextTools::toString(i) + "_";
    newModels.push_back(oldMixture->getNModel(i)->clone());
    newModels.back()->setNamespace(modName);
    allModNames.push_back(modName);
    relprobas.push_back(oldMixture->getNProbability(i));
  }

  for (int i=0; i < nbAdditionalModels; i++) {
  
    NucleotideSubstitutionModel* nucModel = new HKY85(&AlphabetTools::DNA_ALPHABET);
    int currModIdx = oldMixture->getNumberOfModels() + i;
    string modName = "m" + TextTools::toString(currModIdx) + "_";
    allModNames.push_back(modName);
    SubstitutionModel* addedModel = new MutationSelectionModel(standardCode, nucModel, modName);
    newModels.push_back(addedModel);
    
  }

  for (int i=0; i < nbAdditionalModels; i++) {
    relprobas.push_back(60.0/(double)sites_->getNumberOfSites());
  }
  double sum = std::accumulate(relprobas.begin(), relprobas.end(), 0.0);
  for (int i=0; i < relprobas.size(); i++)
    relprobas[i] /= sum;
  vector<double> dummyRates(relprobas.size(), 1.0);
    
  MixtureOfSubstitutionModels * newMixture = new MixtureOfSubstitutionModels(codonAlpha,
									     newModels,
									     relprobas,
									     dummyRates);
  for (int i=0; i < newMixture->getNumberOfModels(); i++) {
    newMixture->getNModel(i)->setNamespace(mixName
					   + TextTools::toString(i + 1)
					   + "_"
					   + allModNames[i]);
  }
  newMixture->setNamespace(mixName);
  
  if (isRoot) {
    newSet->replaceModel(0, newMixture);
  } else {
    newSet->replaceModel(nhtl_->getSubstitutionModelSet()->getModelIndexForNode(bpPosition),
			 newMixture);
  }
  
  vector<string> linkedParamNames({
      "123_HKY85.kappa",
	"123_HKY85.theta",
	"123_HKY85.theta1",
	"123_HKY85.theta2"});

  for (int i=0; i < allModNames.size(); i++) {
    for (string parName : linkedParamNames) {
      string modName = allModNames[i];
      string al = mixName
	+ TextTools::toString(i+1) + "_"
	+ modName
	+ parName + "_"
	  + to_string(newSet->getNumberOfModels());
      string orig = "bp"
	+ TextTools::toString(branchIds_.size())
	+  "_1_m0_"
	+ parName
	+ "_1";
      if (orig != al)
	newSet->aliasParameters(orig, al);
      
    }
    
  }
  
  nhtl_->setSubstitutionModelSet(newSet);
  
}



int BranchMixtureMutSelInference::getCurrentDimension_() {

  int dim = nhtl_->getSubstitutionModelParameters().size()
    + nhtl_->getBranchLengthsParameters().size();

  return dim;

}


double BranchMixtureMutSelInference::aicc_corrected_ (double lnL,
						      size_t dimension,
						      size_t sampleSize) {

  double aicc = (2 * dimension) - (2 * lnL) + (((2 * dimension^2) + (2 * dimension))
					       /
					       (sampleSize - dimension - 1)
					       );

  return aicc;

}


bool BranchMixtureMutSelInference::branchHasBreakpoint_ (int branchId) const {

  if (branchId < 0 || branchId >= tree_->getRootId())
    throw Exception ("BranchMixtureMutSelInference::branchHasBreakpoint_. Invalid branch id requested: " + TextTools::toString(branchId) + ".");
  
  SubstitutionModelSet * oldSet = nhtl_->getSubstitutionModelSet(); 

  vector<int> baseNodes(tree_->getSonsId(tree_->getRootId()));
  bool isBasal = (find(baseNodes.begin(), baseNodes.end(), branchId) != baseNodes.end());
  bool baseHasBp = (oldSet->getModelIndexForNode(baseNodes[0])
		    != oldSet->getModelIndexForNode(baseNodes[1]));
  bool isOldestWithModel =
    (VectorTools::max(oldSet->getNodesWithModel(oldSet->getModelIndexForNode(branchId)))
     == branchId);
  bool nodeHasBp = (isBasal && baseHasBp) || (!isBasal && isOldestWithModel);
  
  return nodeHasBp;
}


double BranchMixtureMutSelInference::likelihoodRatioTest(double lnL1, double lnL0, size_t df) {

  double likelihoodRatio = 2 * (lnL1 - lnL0);

  vector<double> sigVals({0.0001, 0.001, 0.05});

  boost::math::chi_squared dist(df);

  for (int i : sigVals) {
    double upperCritVal = boost::math::quantile(complement(dist, i));
    if (likelihoodRatio > upperCritVal)
      return i;
  }

  return 1.0;

}


void BranchMixtureMutSelInference::newBreakpointAt_ (int newBpNodeIndex) {


  if (tree_->isRoot(newBpNodeIndex))
    throw
      Exception("BranchMixtureMutSelInference::newBreakpointAt_.Breakpoint added at root node.");

  SubstitutionModelSet * oldSet = nhtl_->getSubstitutionModelSet(); 
    
  if  (branchHasBreakpoint_(newBpNodeIndex))
    throw
      Exception("BranchMixtureMutSelinference::newBreakpointAt_.Procedure would overwrite existing breakpoint at " + TextTools::toString(newBpNodeIndex) + ".");

  SubstitutionModelSet * newSet = new SubstitutionModelSet(*oldSet);
  
  StandardGeneticCode * standardCode = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
  CodonAlphabet * codonAlpha = new CodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  NucleotideSubstitutionModel * nucModel = new HKY85(&AlphabetTools::DNA_ALPHABET);
  string mixName = "bp" + TextTools::toString(newBpNodeIndex) + "_";
  string modName = "m0_";
  
  SubstitutionModel * localMutSel = new MutationSelectionModel(standardCode, nucModel, modName); 
  
  vector<SubstitutionModel*> startModel;
  startModel.push_back(localMutSel);
  MixtureOfSubstitutionModels * localMixture = new MixtureOfSubstitutionModels(codonAlpha, startModel);
  
  localMixture->setNamespace(mixName);
  for (int i=0; i < localMixture->getNumberOfModels(); i++) {
    localMixture->getNModel(i)->setNamespace(mixName + TextTools::toString(i+1) + "_" + modName);
  }

  vector <int> newBpBranches = getDescendantsWithModel(*newSet, newBpNodeIndex); 
      
  newSet->addModel(localMixture, newBpBranches);

  nhtl_->setSubstitutionModelSet(newSet);
								    
}


vector<int> BranchMixtureMutSelInference::getDescendantsWithModel(const SubstitutionModelSet& mset,
								  int nodeId) {

  vector <int> descendants;
  getDescendantsWithModel_(mset, nodeId, &descendants);
  return descendants;
  
}


void BranchMixtureMutSelInference::getDescendantsWithModel_(const SubstitutionModelSet& mset,
							    int nodeId,
							    vector<int>* descendants) {
  descendants->push_back(nodeId);

  int modelIndex = mset.getModelIndexForNode(nodeId);
  vector <int> nodesWithModel = mset.getNodesWithModel(modelIndex);

  vector <int> sons = tree_->getSonsId(nodeId);
  for (auto sonId = sons.begin(); sonId != sons.end(); ++sonId) {
    vector <int>::iterator hasModel = find(nodesWithModel.begin(), nodesWithModel.end(), *sonId);
    if (hasModel != nodesWithModel.end()) {
      getDescendantsWithModel_(mset, *sonId, descendants);
    }
  }

}

// Unused at present.
ParameterList BranchMixtureMutSelInference::makeBackgroundParameterList(const VVdouble& bgParams,
									const Vdouble& bgWeights) {
 

  ParameterList initialParams;

  vector<int> bps;
  bps.push_back(tree_->getRootId());
  for (int i=0; i < assumedBreakpoints_.size(); i++)
    bps.push_back(assumedBreakpoints_[i]);

  for (int bp=0; bp < bps.size(); bp++) {
  
    for (int i=0; i < bgParams.size(); i++) {
      Parameter relprob("bp" + TextTools::toString(bps[bp])
			+ "_relproba" + TextTools::toString(i) + "_"
			+ TextTools::toString(bp + 1),
			bgWeights[i]);
      initialParams.addParameter(relprob);
      
      for (int j=0; j < bgParams[i].size(); j++) {
	Parameter p("bp" + TextTools::toString(bps[bp]) + "_"
		    + TextTools::toString(i + 1)
		    + "_m" + TextTools::toString(i)
		    + "_fit" + TextTools::toString(j) + "_"
		    + TextTools::toString(bp + 1),
		    bgParams[i][j]);
	initialParams.addParameter(p);
      }
      
    }
    
  }

  return initialParams;

}


/**********************************************************************************************/
