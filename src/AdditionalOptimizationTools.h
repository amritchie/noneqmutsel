// File: AdditionalOptimizationTools.h
//
//

/**
 *Extra optimization procedures required
 * for non-equilibrium mutation-selection
 * analyses.
 *
 */

#ifndef _ADDITIONAL_OPTIMIZATIONTOOLS_H_
#define _ADDITIONAL_OPTIMIZATIONTOOLS_H_

#include <Bpp/Phyl/Likelihood/NonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Numeric/Function/Optimizer.h>

#include "MixedSubstitutionModelFunction.h"

namespace bpp
{
  
  class AdditionalOptimizationTools {

  public:
    AdditionalOptimizationTools();
    virtual ~AdditionalOptimizationTools();

  public:
    static unsigned int optimizeNumericalParametersNHMixed
      (NonHomogeneousTreeLikelihood* nhtl,
       const ParameterList& parameters,
       OptimizationListener* listener = 0,
       unsigned int nstep = 1,
       double tolerance = 0.000001,
       unsigned int tlEvalMax = 1000000,
       OutputStream* messageHandler = ApplicationTools::message.get(),
       OutputStream* profiler = ApplicationTools::message.get(),
       bool reparametrization = false,
       unsigned int verbose = 1,
       const std::string& optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON,
       const std::string& optMethodModel = OptimizationTools::OPTIMIZATION_BRENT);


    static unsigned int optimizeMixedSubstitutionModelFunction(MixedSubstitutionModelFunction* mixfunc,   
						      const ParameterList& parameters,
						      OptimizationListener* listener,
						      unsigned int nstep,
						      double tolerance,
						      unsigned int tlEvalMax,
						      OutputStream* messageHandler,
						      OutputStream* profiler,
						      unsigned int verbose);
    
  };


} // namespace bpp

#endif //_ADDITIONAL_OPTIMIZATIONTOOLS_H_
