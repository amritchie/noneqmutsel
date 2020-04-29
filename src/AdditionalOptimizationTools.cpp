//
// File: AdditionalOptimizationTools.cpp
//
//

// This is really a scratchpad for messing around with Bpp's optimisers
// The optimisation actually works fine, though it could be faster; it's the sample size
// that is causing all the problems. Probably no further need to mess with this.

#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>

//#include "/home/aritchie/src/nemss2/src/RNonHomogeneousMultiProfileTreeLikelihood.h"
#include "MixedSubstitutionModelFunction.h"
#include "AdditionalOptimizationTools.h"

using namespace bpp;

using namespace std;


AdditionalOptimizationTools::AdditionalOptimizationTools() {}

AdditionalOptimizationTools::~AdditionalOptimizationTools() {}


unsigned int AdditionalOptimizationTools::optimizeNumericalParametersNHMixed(NonHomogeneousTreeLikelihood* nhtl,
								      const ParameterList& parameters,
								      OptimizationListener* listener,
								      unsigned int nstep,
								      double tolerance,
								      unsigned int tlEvalMax,
								      OutputStream* messageHandler,
								      OutputStream* profiler,
								      bool reparametrization,
								      unsigned int verbose,
								      const std::string& optMethodDeriv,
								      const std::string& optMethodModel)
{
  DerivableSecondOrder* f = nhtl;
  ParameterList pl = parameters;

  // Shall we reparametrize the function to remove constraints?
  unique_ptr<DerivableSecondOrder> frep;
  if (reparametrization)
    {
      frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, parameters));
      f = frep.get();

      // Reset parameters to remove constraints:
      pl = f->getParameters().subList(parameters.getParameterNames());
    }

  // ///////////////
  // Build optimizer:

  // Branch lengths

  MetaOptimizerInfos* desc = new MetaOptimizerInfos();
  MetaOptimizer* poptimizer = 0;
  AbstractNumericalDerivative* fnum = new ThreePointsNumericalDerivative(f);
  fnum->setInterval(0.000001);

   if (optMethodDeriv == OptimizationTools::OPTIMIZATION_GRADIENT)
     desc->addOptimizer("Branch length parameters", new ConjugateGradientMultiDimensions(f), nhtl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
   else if (optMethodDeriv == OptimizationTools::OPTIMIZATION_NEWTON)
     desc->addOptimizer("Branch length parameters", new PseudoNewtonOptimizer(f), nhtl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
   else if (optMethodDeriv == OptimizationTools::OPTIMIZATION_BFGS)
     desc->addOptimizer("Branch length parameters", new BfgsMultiDimensions(f), nhtl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
   else
     throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);

  // Other parameters

  if (optMethodModel == OptimizationTools::OPTIMIZATION_BRENT)
    {
      ParameterList plsm = parameters.getCommonParametersWith(nhtl->getSubstitutionModelParameters());
      desc->addOptimizer("Substitution model parameter", new SimpleMultiDimensions(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);


      ParameterList plrd = parameters.getCommonParametersWith(nhtl->getRootFrequenciesParameters());
      desc->addOptimizer("Profile frequency parameters", new SimpleMultiDimensions(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
      poptimizer = new MetaOptimizer(f, desc, nstep);
    }
  else if (optMethodModel == OptimizationTools::OPTIMIZATION_BFGS)
    {
      vector<string> vNameDer;

      ParameterList plsm = parameters.getCommonParametersWith(nhtl->getSubstitutionModelParameters());
      ParameterList plnm;
      for (int p = 0; p < plsm.size(); p++) {
	int found = plsm[p].getName().find("theta");
	if (found != string::npos) {
	  plnm.addParameter(plsm[p]);
	}
      }

      for (int p = 0; p < plnm.size(); p++) {
	plsm.deleteParameter(plnm[p].getName());
      }

      vNameDer = plsm.getParameterNames();

      //ParameterList plrd = parameters.getCommonParametersWith(nhtl->getProfileFrequencyParameters());

      //vector<string> vNameDer2 = plrd.getParameterNames();

      //vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
      fnum->setParametersToDerivate(vNameDer);

      
      desc->addOptimizer("Nucleotide model parameters", new SimpleMultiDimensions(f), plnm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
      
      //desc->addOptimizer("Profile frequency parameters", new SimpleMultiDimensions(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_FULL);
      
      desc->addOptimizer("stuffup", new DownhillSimplexMethod(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_FULL);
      poptimizer = new MetaOptimizer(f, desc, nstep);
      
      // desc->addOptimizer("Model & profile frequency parameters", new BfgsMultiDimensions(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
      //poptimizer = new MetaOptimizer(fnum, desc, nstep);
    }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodModel);

  poptimizer->setVerbose(verbose);
  poptimizer->setProfiler(profiler);
  poptimizer->setMessageHandler(messageHandler);
  poptimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  poptimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  NaNListener* nanListener = new NaNListener(poptimizer, nhtl);
  poptimizer->addOptimizationListener(nanListener);
  if (listener)
    poptimizer->addOptimizationListener(listener);
  poptimizer->init(pl);
  poptimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  unsigned int nb = poptimizer->getNumberOfEvaluations();
  delete poptimizer;
  return nb;
}


// Don't use this; it was a dead-end attempt to make the mixture optimisation faster.
unsigned int AdditionalOptimizationTools::optimizeMixedSubstitutionModelFunction(MixedSubstitutionModelFunction* mixfunc,   
							    const ParameterList& parameters,
							    OptimizationListener* listener,
							    unsigned int nstep,
							    double tolerance,
						            unsigned int tlEvalMax,
							    OutputStream* messageHandler,
							    OutputStream* profiler,
							    unsigned int verbose){
  Function * f = mixfunc;
  ParameterList pl = parameters;

  
  MetaOptimizerInfos* desc = new MetaOptimizerInfos();
  MetaOptimizer* poptimizer = 0;
  
  // Branch lengths

  // desc->addOptimizer("Branch length parameters",
  // 		     new PseudoNewtonOptimizer(f),
  // 		     mixfunc->getBranchLengthsParameters().getParameterNames(),
  // 		     2,
  // 		     MetaOptimizerInfos::IT_TYPE_FULL);

  // Frequency parameters


  
  desc->addOptimizer("Frequency parameters",
		     new SimpleMultiDimensions(f),
		     parameters.getParameterNames(),
		     0,
		     MetaOptimizerInfos::IT_TYPE_STEP);

  poptimizer = new MetaOptimizer(f, desc, nstep);

  poptimizer->setVerbose(verbose);
  poptimizer->setProfiler(profiler);
  poptimizer->setMessageHandler(messageHandler);
  poptimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  poptimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  NaNListener* nanListener = new NaNListener(poptimizer, f);
  poptimizer->addOptimizationListener(nanListener);
  if (listener)
    poptimizer->addOptimizationListener(listener);
  poptimizer->init(pl);
  poptimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  unsigned int nb = poptimizer->getNumberOfEvaluations();
  delete poptimizer;
  return nb;

  
}
