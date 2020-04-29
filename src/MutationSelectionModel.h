/**File: MutationSelectionModel.h
 **
 **Substitution model class for simulating under a time-reversible Mutation-Selection Model
 **(see Halpern & Bruno, 1998).
 **/

#ifndef _MUTATION_SELECTION_MODEL_H_
#define _MUTATION_SELECTION_MODEL_H_

#include <Bpp/Phyl/Model/Codon/AbstractCodonSubstitutionModel.h>

namespace bpp
{

  class MutationSelectionModel:
    public AbstractCodonSubstitutionModel
  {
  private:

    double fixprob_[20][20];
    double scaleFactor_;
    
  protected:

    std::map<int, double> fitset_;
    double popsize_;
    StateMap* stateMap_;
    bool neutrality_;
    bool logTransformedFitnesses_;
    
    void updateFixationProbabilities_();
    void normalizeFitnesses_();
    double updateScale_() const;
    
  public:
  
    MutationSelectionModel(const GeneticCode* gCode,
			   NucleotideSubstitutionModel* pmod,
			   const std::string& prefix,
			   bool logTransformedFitnesses=false,
			   bool neutrality=false);

    MutationSelectionModel(const GeneticCode* gCode,
			   NucleotideSubstitutionModel* pmod,
			   const std::map<int, double>& fitset,
			   const std::string& prefix,
			   bool logTransformedFitnesses=false,
			   bool neutrality=false);

  MutationSelectionModel(const MutationSelectionModel& mutsel) :
    AbstractParameterAliasable::AbstractParameterAliasable(mutsel),
      AbstractCodonSubstitutionModel::AbstractCodonSubstitutionModel(mutsel),
      fitset_(mutsel.fitset_),
      popsize_(mutsel.popsize_),
      stateMap_(mutsel.stateMap_),
      logTransformedFitnesses_(mutsel.logTransformedFitnesses_),
      neutrality_(mutsel.neutrality_)
	{}
    
    MutationSelectionModel* clone() const {
      return new MutationSelectionModel(*this);
    }

    virtual void fireParameterChanged(const ParameterList& parameters);

    double getCodonsMulRate(size_t i, size_t j) const;
  
    const FrequenciesSet* getFrequenciesSet() const;

    std::string getName() const {
      return "MutSel";
    }

    double getScale() const;

    double computeFixationProbability(double fit_source, double fit_target);

  };

} //namespace bpp

#endif //_MUTATION_SELECTION_MODEL_H_
