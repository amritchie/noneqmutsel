

/**File: MutationSelectionModel.cpp
 **
 **Substitution model class for simulating under a nonreversible Mutation-Selection Model
 **(see Halpern & Bruno, 1998).
 **
 **/

#include <math.h>
#include <algorithm>
#include <iterator>

#include <Bpp/Phyl/Model/Codon/AbstractCodonSubstitutionModel.h>
#include <Bpp/Phyl/Model/StateMap.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Utils/MapTools.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Numeric/Parameter.h>

#include "MutationSelectionModel.h"

using namespace bpp;
using namespace std;

MutationSelectionModel::MutationSelectionModel(const GeneticCode* gCode,
					       NucleotideSubstitutionModel* pmod,
					       const string& prefix,
					       bool logTransformedFitnesses,
					       bool neutrality
					       ) :
  AbstractCodonSubstitutionModel(gCode, pmod, prefix, false),
  AbstractParameterAliasable(prefix),
  stateMap_(new CanonicalStateMap(gCode->getSourceAlphabet(), false))
{
  popsize_ = 10000.0;
  neutrality_ = neutrality;
  logTransformedFitnesses_ = logTransformedFitnesses;
  vector<double> fitvec(20);
  double initVal = logTransformedFitnesses? 0.0 : 1.0;
  if (neutrality_) {
    VectorTools::fill(fitvec, initVal);
  } else {
    VectorTools::fill(fitvec, initVal);
    //vector<double> initvec({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    //for (int i=0;i<initvec.size();i++) fitvec[i] = initvec[i];
  }
  for (int i=0; i < fitvec.size(); i++) {
    fitset_.insert(pair<int, double>(i, fitvec[i]));
  }
  for(int i=1; i<fitset_.size(); i++) {
    string parName = getNamespace() + "fit" + to_string(i);
    if (logTransformedFitnesses) {
      addParameter_(new Parameter(parName, fitset_[i]));
    } else {
      addParameter_(new Parameter(parName, fitset_[i], &Parameter::R_PLUS, false));
    }
  }
  updateFixationProbabilities_();
  computeFrequencies(true);
  //if (!neutrality_) updateScale_();
  updateMatrices();
}

MutationSelectionModel::MutationSelectionModel(const GeneticCode* gCode,
					       NucleotideSubstitutionModel* pmod,
					       const map<int, double>& fitset,
					       const string& prefix,
					       bool logTransformedFitnesses,
					       bool neutrality) :
  AbstractCodonSubstitutionModel(gCode, pmod, prefix, false),
  AbstractParameterAliasable(prefix)
{
  popsize_ = 10000;
  neutrality_= neutrality;
  logTransformedFitnesses_ = logTransformedFitnesses;
  size_t nStates = gCode->getTargetAlphabet()->getNumberOfStates();
  if (fitset.size() != (nStates-1)) {
    throw DimensionException("Fitness map incorrectly sized", fitset.size(), (nStates-1));
  }
  copy(fitset.begin(), fitset.end(), inserter(fitset_, fitset_.begin()));
  for (int i=0; i < nStates; i++) {
    if (fitset_.find(i) == fitset_.end()) {
      fitset_.insert(pair<int, double>(i,0));
      break;
    }
  }
  //normalizeFitnesses_();
  updateFixationProbabilities_();
  for(auto p : fitset_) {
    string parName = getNamespace() + "fit" + to_string(p.first);
    Parameter prm(parName, p.second, NULL, false);
    addParameter_(&prm);
  }
  //if (!neutrality_) updateScale_();
}

void MutationSelectionModel::fireParameterChanged(const ParameterList& parameters) {
  for (int i=1; i < fitset_.size(); i++) {
    string parName = "fit" + to_string(i);
    if (hasParameter(parName)) {
      fitset_[i] = getParameter(parName).getValue();
    }
  }
  updateFixationProbabilities_();

  // bool rescale = false;

  // for (int i=0; i < parameters.size(); i++) {

  //   size_t isfreq = parameters[i].getName().find("theta");
  //   size_t israte = parameters[i].getName().find("kappa");
  //   if (isfreq != string::npos || israte != string::npos) {
  //     rescale = true;
  //     break;
  //   }
  // }

  // if (rescale) scaleFactor_ = updateScale_();
  
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double MutationSelectionModel::getCodonsMulRate (size_t i, size_t j) const {
  const GeneticCode* code = getGeneticCode();
  int aa_i = code->translate(stateMap_->getAlphabetStateAsInt(i));
  int aa_j = code->translate(stateMap_->getAlphabetStateAsInt(j));                                  
  double fixProb = fixprob_[aa_i][aa_j]; 
  return fixProb;
}

void MutationSelectionModel::updateFixationProbabilities_() {
  for (int i=0; i < 20; i++) {
    for (int j=0; j < 20; j++) {
      if (i == j) {
	fixprob_[i][j] = 1.0;
      } else {
	fixprob_[i][j] = (2.0 * popsize_) * computeFixationProbability(fitset_[i], fitset_[j]);	
      }
    }
  }
}

double MutationSelectionModel::computeFixationProbability(double fit_source, double fit_target) {
  double s = logTransformedFitnesses_? (fit_target - fit_source) : (fit_target/fit_source - 1.0);
  double fab = logTransformedFitnesses_?
    ((1.0 - exp(-1.0 * s * (1.0/popsize_)))/(1.0 - exp(-2.0 * s))):
    ((1.0 - exp(-1.0 * s)) / (1.0 - exp(-2.0 * popsize_ * s))); //s/(1.0 - exp(-1.0*s));
  if (std::fabs(s) < NumConstants::TINY()) {
    return 1.0/popsize_;
  } else {
    // if (std::fabs(fab) < NumConstants::VERY_TINY()) {
    //   cout << fit_source << " " << fit_target << endl;
    //   throw
    // 	Exception("WE UNDERFLEW!!!");
    // }
    return fab;
  }
}

const FrequenciesSet* MutationSelectionModel::getFrequenciesSet() const {
  return 0;
}

double MutationSelectionModel::updateScale_() const {
  if (neutrality_) {
    return 1.0;
  } else {
    SubstitutionModel* nmod = VSubMod_[0]->clone();
    MutationSelectionModel* scalingModel
      = new MutationSelectionModel(getGeneticCode(),
				   dynamic_cast<NucleotideSubstitutionModel*>(nmod),
				   getNamespace(),
				   logTransformedFitnesses_,
				   true);

    RowMatrix<double> scalingMatrix = scalingModel->getGenerator();
    vector<double> neutralFreqs = scalingModel->getFrequencies();
    vector<double> v;
    MatrixTools::diag(scalingMatrix, v);
    double scale =  -VectorTools::scalar<double, double>(v, neutralFreqs);
    
    delete scalingModel;

    return scale;		   
  }
}

double MutationSelectionModel::getScale() const {
  return updateScale_();
}
