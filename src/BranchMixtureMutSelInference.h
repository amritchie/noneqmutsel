
/** File : BranchMixtureMutSelInference.h **/

#ifndef _BRANCH_MIXTURE_MUT_SEL_INFERENCE_H_
#define _BRANCH_MIXTURE_MUT_SEL_INFERENCE_H_

#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>

namespace bpp
{

  class BranchMixtureMutSelInference {


  protected:

    const Tree * tree_;
    const SiteContainer * sites_;

    size_t maxBreakpoints_;
    size_t sampleSize_;

    std::vector<int> assumedBreakpoints_;
    std::vector<int> inferredBreakpoints_;
    std::vector<int> branchIds_;

    std::map<int, vector<double> > positionLikelihoods_;
    std::map<int, double> aiccs_;
    std::vector < std::pair <string, std::time_t> > timings_;
    
    DRNonHomogeneousTreeLikelihood * nhtl_;
    DRNonHomogeneousTreeLikelihood * storedLikelihood_;

    
  public:

    BranchMixtureMutSelInference (const Tree * tree,
				  const SiteContainer * sites,
				  const std::vector<int>& assumedBreakpoints);

    BranchMixtureMutSelInference (const BranchMixtureMutSelInference& inf);

    BranchMixtureMutSelInference& operator= (const BranchMixtureMutSelInference& inf);

    virtual ~BranchMixtureMutSelInference ();

    /**********************************************************************************************/
    
    void bppInferBackgroundMixture(size_t maxInitialModels,
				   bool nsRootFreqs = false);

    void bppFindSiteFrequencyShifts(const ParameterList& initialParams,
				    size_t nBackgroundCats,
				    size_t maxFreqShifts);

    double bppTestSpecificBranch(int testBranchIndex,
				 size_t MaxInitialModels,
				 size_t maxBreakpointModels);

    double bppTestBestBranch(size_t maxBreakpoints,
			     size_t maxInitialModels,
			     size_t maxBreakpointModels);

    void writeParameterEstimates(string path);

    void writeBreakpointEstimates(string path);

    void writeEquilibriumFrequencies(string path);

    void writeTimings(string path);

    void setAssumedBreakpoints(const vector<int>& breakpoints);

    void getAllMixtureComponents(vector <SubstitutionModel*>& models) const;

    double computeLogLikelihood();

    double computeLogLikelihoodForParameters(const ParameterList& pl, int nbCategories);

    Vdouble computeSiteLogLikelihoodsForParameters(const ParameterList & pl, int nbCategories);

    void getEmpiricalCodonFreqs(VVdouble& tipFreqs, Vdouble& averageFreqs);

    void computeTipDistrsForParameters(const ParameterList& pl,
				       int nbCategories,
				       VVdouble& outDists);
    
    /**********************************************************************************************/
    
  private:

    // Likelihood initialisation / reversion
    
    void init_();

    void storeLikelihood_();

    void revertLikelihood_();
    
    void resetLikelihood_();

    // Higher level optimization procedures

    int optimizeBreakpointPosition_(size_t nbAddedModels, bool optimizeParameters);
    
    double optimizeBPMixture_(int bpPosition, size_t maxModels, bool optimizeParameters);

    double optimizeMixtureParams_(int bpPosition,
				  double tolerance,
				  vector <size_t> submodels,
				  bool includeBranchLengths);

    double optimizeAllParameters_();

    double optimizeParams_(ParameterList optParams, double tolerance);

    double optimizeSiteFrequenciesAllBreakpoints_(); 


    // Model set modification

    void addModelsToBp_(int bpPosition,
			size_t nbAdditionalModels);
    

    // Utility procedures

    ParameterList * getParametersForBPSubmodels_(int bpPosition,
						 vector<size_t> submodels);

    ParameterList * getFrequencyParamsForBP_(int bpPosition);

    int getCurrentDimension_();
    
    double aicc_corrected_ (double lnL, size_t dimension, size_t sampleSize);

    double likelihoodRatioTest(double lnL1, double lnL2, size_t df);
    
    bool branchHasBreakpoint_(int branchId) const;
    
    void computeTipDistrs_(int branch, Vdouble * startDist, VVdouble& outDists);

    void newBreakpointAt_(int newBpNodeIndex);

    vector <int> getDescendantsWithModel(const SubstitutionModelSet& mset,
					 int nodeId);

    void getDescendantsWithModel_(const SubstitutionModelSet& mset,
				  int nodeId,
				  vector <int> * descendants);
    
    ParameterList makeBackgroundParameterList(const VVdouble& bgParams,
					      const Vdouble& bgWeights);
    
  };



} //namespace bpp;

#endif // _BRANCH_MIXTURE_MUT_SEL_INFERENCE_H_
