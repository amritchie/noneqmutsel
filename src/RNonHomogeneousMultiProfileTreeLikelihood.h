/** File: RNonHomogeneousMultiProfileTreeLikelihood.h **/

#ifndef _R_NON_HOMOGENEOUS_MULTI_PROFILE_TREE_LIKELIHOOD_H_
#define _R_NON_HOMOGENEOUS_MULTI_PROFILE_TREE_LIKELIHOOD_H_

#include <Bpp/Phyl/Likelihood/AbstractTreeLikelihood.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>

//#include "/home/aritchie/src/nemss2/src/NonHomogeneousSiteProfileMixture.h"

namespace bpp
{
  class RNonHomogeneousMultiProfileTreeLikelihood:
    public AbstractTreeLikelihood
  {

  public:

  /* !TODO: There is probably no sensible way for these iterators to work
   * for profile mixtures like this.
   * but only the substitution mapping methods use them, so this will do for now.
   * 
   */

    class ConstNoPartitionMixedBranchModelIterator:
    public ConstBranchModelIterator
    {
    private:
      std::vector<ConstNoPartitionBranchModelDescription> branchModelDescriptions_;
      size_t index_;
      size_t nbModels_;
      size_t nbSites_;
     
    public:

    ConstNoPartitionMixedBranchModelIterator(const std::vector <SubstitutionModelSet*>& profileSet,
					     size_t nodeId,
					     size_t nbSites):
      branchModelDescriptions_(),
	index_(0),
	nbModels_(profileSet.size()),
	nbSites_(nbSites)
	  {
	    for (size_t i=0; i<nbModels_; ++i) {
	      branchModelDescriptions_.push_back(ConstNoPartitionBranchModelDescription(profileSet[i]->getModelForNode(nodeId), nbSites_));
	    }
	  }

      ConstNoPartitionBranchModelDescription* next() {
	if (!hasNext())
	  throw Exception("AbstractNonHomogeneousTreeLikelihood::ConstNoPartitionMixedBranchModelIterator::next(). No more branch models in the set.");
	return &branchModelDescriptions_[index_++];
      }

      bool hasNext() const { return index_ < nbModels_; }
     
    };
 
  class ConstNonHomogeneousMixedSiteModelIterator:
    public ConstSiteModelIterator
  {
  private:
    std::vector<ConstNoPartitionSiteModelDescription> siteModelDescriptions_;
    size_t index_;
    size_t nbProfiles_;
    size_t modelsPerProfile_;
    size_t nbModels_;

  public:
  ConstNonHomogeneousMixedSiteModelIterator(const std::vector<SubstitutionModelSet*>& profileSet) :
    siteModelDescriptions_(),
      index_(0),
      nbProfiles_(profileSet.size()),
      modelsPerProfile_(profileSet[0]->getNumberOfModels()),
      nbModels_(0)
	{
	  nbModels_ = nbProfiles_ * modelsPerProfile_;
	  for (size_t i = 0; i < nbProfiles_; ++i) {
	    for (size_t j = 0; j < modelsPerProfile_; ++j) {
	      siteModelDescriptions_.push_back(ConstNoPartitionSiteModelDescription(profileSet[i]->getModel(j), profileSet[i]->getNodesWithModel(j)));
	    }
	  }
	}

  public:
    ConstSiteModelDescription* next()
    {
      if (!hasNext())
	throw Exception("RNonHomogeneousMixedProfileTreeLikelihood::ConstNonHomogeneousMixedSiteModelIterator::next(). No more sites in the set.");
      return &siteModelDescriptions_[index_++];
    }

    bool hasNext() const { return index_ < nbModels_; }
  };

 
 private:

  double minusLogLik_;
   
  mutable DRASRTreeLikelihoodData* likelihoodData_;
  

 protected:

  std::vector <SubstitutionModelSet*>  profileSet_;

  SimpleDiscreteDistribution* profileFreqs_;
  
  ParameterList brLenParameters_;

  mutable std::map<int, VVVdouble> pxy_;

  mutable std::map<int, VVVdouble> dpxy_;

  mutable std::map<int, VVVdouble> d2pxy_;

  VVdouble rootFreqs_;

  std::vector<Node*> nodes_;

  mutable std::map<int, const Node*> idToNode_;

  size_t nbSites_,
    nbDistinctSites_,
    nbProfiles_,
    nbStates_,
    nbNodes_;

  double minimumBrLen_;
  double maximumBrLen_;
  std::unique_ptr<Constraint> brLenConstraint_;

  int root1_, root2_;

 public:

  RNonHomogeneousMultiProfileTreeLikelihood(const Tree& tree,
					    const SiteContainer& sites,
					    std::vector<SubstitutionModelSet*>& profileSet,
					    SimpleDiscreteDistribution* profileFreqs);

  
    /**
     * @brief Copy constructor
     *
     * This constructor is to be called by the derived class copy constructor.
     */

  RNonHomogeneousMultiProfileTreeLikelihood(const RNonHomogeneousMultiProfileTreeLikelihood& lik);

  
  /**
   * @brief Assignment operator
   *
   * This operator is to be called by the derived class operator.
   */

  RNonHomogeneousMultiProfileTreeLikelihood& operator=(const RNonHomogeneousMultiProfileTreeLikelihood& lik);

  virtual ~RNonHomogeneousMultiProfileTreeLikelihood() {}

  RNonHomogeneousMultiProfileTreeLikelihood* clone() const { return new RNonHomogeneousMultiProfileTreeLikelihood(*this); }

 private:

  void init_ (const Tree& tree,
	      std::vector <SubstitutionModelSet*>& modelSet,
	      SimpleDiscreteDistribution* rDist);
  
 public:

  /**
   * TreeLikelihood interface members
   *
   *
   */

  size_t getNumberOfStates() const { return profileSet_[0]->getNumberOfStates(); } 
    
  const std::vector<int>& getAlphabetStates() const { return profileSet_[0]->getAlphabetStates(); } 
    
  int getAlphabetStateAsInt(size_t i) const { return profileSet_[0]->getAlphabetStateAsInt(i); }
  
  std::string getAlphabetStateAsChar(size_t i) const { return profileSet_[0]->getAlphabetStateAsChar(i); }
  
  void initialize();
    
  ParameterList getBranchLengthsParameters() const;
    
  ParameterList getSubstitutionModelParameters() const;
  
  ParameterList getProfileFrequencyParameters() const;
  
  std::vector<TransitionModel*> getModelsForNode(int nodeId) const;

  const VVdouble& getRootFrequenciesPerProfile() const {return rootFreqs_;}
  
  const std::vector<double>& getRootFrequencies(size_t siteIndex) const;
  
  ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const { return new ConstNoPartitionMixedBranchModelIterator(profileSet_, static_cast<size_t>(nodeId), nbDistinctSites_);    
  }
  
  ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const {
    return new ConstNonHomogeneousMixedSiteModelIterator(profileSet_);
  }

  double getLikelihoodForASiteForAState (size_t site, int state) const;
  
  double getLogLikelihoodForASiteForAState (size_t site, int state) const;
  
  TransitionModel* getModelForSite(int nodeId, size_t siteIndex);
  
  const TransitionModel* getModelForSite(int nodeId, size_t siteIndex) const;

  VVVdouble getTransitionProbabilitiesPerProfile(int nodeId, size_t siteIndex) const { return pxy_[nodeId]; }
  
  VVdouble getTransitionProbabilities(int nodeId, size_t siteIndex) const;

  ParameterList getDerivableParameters() const;
  
  ParameterList getNonDerivableParameters() const;
  
  
  /**
   *Replacements for methods implemented in RNonHomogeneousTreeLikelihood
   *
   */

  void setData(const SiteContainer& sites);
  double getLikelihood() const;
  double getLogLikelihood() const;
  double getLikelihoodForASite(size_t site) const;
  double getLogLikelihoodForASite(size_t site) const;
  size_t getSiteIndex(size_t site) const { return likelihoodData_->getRootArrayPosition(site); }

  double getLikelihoodForASiteForAProfile(size_t site, size_t profile) const;
  double getLogLikelihoodForASiteForAProfile(size_t site, size_t profile) const;
  double getLikelihoodForASiteForAProfileForAState(size_t site, size_t profile, int state) const;
  double getLogLikelihoodForASiteForAProfileForAState(size_t site, size_t profile, int state) const;

  DRASRTreeLikelihoodData* getLikelihoodData() {return likelihoodData_; } 

  const DRASRTreeLikelihoodData* getLikelihoodData() const {return likelihoodData_; }

  virtual void computeTreeLikelihood();

  virtual double getDLikelihoodForASiteForAProfile(size_t site, size_t profile) const;

  virtual double getDLikelihoodForASite(size_t site) const;

  virtual double getDLogLikelihoodForASite(size_t site) const;
		
  virtual double getDLogLikelihood() const;
		
  virtual void computeTreeDLikelihood(const std::string& variable);

  virtual double getD2LikelihoodForASiteForAProfile(size_t site, size_t rateClass) const;

  virtual double getD2LikelihoodForASite(size_t site) const;

  virtual double getD2LogLikelihoodForASite(size_t site) const;
		
  virtual double getD2LogLikelihood() const;
		
  virtual void computeTreeD2Likelihood(const std::string& variable);
  
  /**
   * Replacements for the NonHomogeneousTreeLikelihood members.
   *
   */

  const std::vector <SubstitutionModelSet*>&  getProfileSet() const {return profileSet_;}

  SimpleDiscreteDistribution* getProfileFreqs() const {return profileFreqs_;}
  
  void setProfiles(std::vector <SubstitutionModelSet*>& modelSet);


  /**
   * Replacements for the AbstractNonHomogeneousTreeLikelihood members.
   *
   */

  virtual void initParameters();
  virtual void applyParameters();
  virtual void initBranchLengthsParameters();

  virtual void setMinimumBranchLength(double minimum);
  virtual void setMaximumBranchLength(double maximum);

  virtual double getMinimumBranchLength() const { return minimumBrLen_; }
  virtual double getMaximumBranchLength() const { return maximumBrLen_; }
  


  /**
   * Function interface members.
   *
   */

  void setParameters(const ParameterList& parameters);
  double getValue() const;

  double getFirstOrderDerivative(const std::string& variable) const;
  
  virtual double getSecondOrderDerivative (const std::string &variable) const;

  virtual double getSecondOrderDerivative (const std::string &variable1,
					   const std::string &variable2) const {return 0;};



 protected:

  virtual void computeTreeLikelihoodForProfiles(std::vector <size_t>& dirtyProfiles);
  
  virtual void computeSubtreeLikelihood(const Node* node);

  virtual void computeSubtreeLikelihoodForProfiles(const Node* node,
						    std::vector <size_t>& dirtyProfiles);

  virtual void computeTransitionProbabilitiesForNodeForProfile(const Node* node, size_t profile);
  
  virtual void computeDownSubtreeDLikelihood(const Node* node);

  virtual void computeDownSubtreeD2Likelihood(const Node* node);

  void fireParameterChanged(const ParameterList & params);


 protected:

  virtual void computeAllTransitionProbabilities ();

  virtual void computeTransitionProbabilitiesForNode (const Node*);

};

}
#endif //_R_NON_HOMOGENEOUS_MULTI_PROFILE_TREE_LIKELIHOOD_H_
