/** File: NonHomogeneousSiteProfileMixture.h **/

#ifndef _NON_HOMOGENEOUS_SITE_PROFILE_MIXTURE_H_
#define _NON_HOMOGENEOUS_SITE_PROFILE_MIXTURE_H_

#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

namespace bpp
{

  class  NonHomogeneousSiteProfileMixture :
    public AbstractParameterAliasable
  {

  protected:

    std::vector <SubstitutionModelSet *> profileSet_;

    std::vector <ParameterList> profileParameters_;
    
    
  public:

    /* Constructors, destructors and operator overloads */
    
    NonHomogeneousSiteProfileMixture (std::vector <SubstitutionModelSet*>& profileSet);

    NonHomogeneousSiteProfileMixture (const NonHomogeneousSiteProfileMixture& mix);

    NonHomogeneousSiteProfileMixture& operator=(const NonHomogeneousSiteProfileMixture& mix);

    NonHomogeneousSiteProfileMixture* clone() const {
      return new NonHomogeneousSiteProfileMixture(*this);
    }
    
    virtual ~NonHomogeneousSiteProfileMixture () {}
    

  public:

    /* Methods related to the profile vector */

    const SubstitutionModelSet * getProfile (size_t k) const { return profileSet_[k]; }

    int getNumberOfProfiles () const { return profileSet_.size(); }

    void addProfile(SubstitutionModelSet *  profile);

    void aliasParameterAcrossProfiles(std::string name);

    bool allAreFullySetUpFor(const Tree& tree) const;
      
    
    /* Utility methods carried over from SubstitutionModelSet... very ugly! */
    
    std::vector <int> getNodesWithParameter(std::string pname) const;

    const TransitionModel * getModelForProfileForNode(size_t profile, int nodeId) const;

    const SubstitutionModel * getSubstitutionModelForProfileForNode(size_t profile, int nodeId) const;
    
    TransitionModel * getModelForProfileForNode(size_t profile, int nodeId);

    SubstitutionModel * getSubstitutionModelForProfileForNode(size_t profile, int nodeId);

    std::vector <TransitionModel * > getAllModelsForNode(int nodeId);

    std::vector <SubstitutionModel * > getAllSubstitutionModelsForNode(int nodeId);

    ParameterList getNodeParametersForProfile(size_t profile) const;

    const std::vector <double> getRootFrequenciesForProfile(size_t profile) const;

    size_t getNumberOfStates() const {return profileSet_[0]->getNumberOfStates();}

    const Alphabet * getAlphabet() const {return profileSet_[0]->getAlphabet(); }
    
    const std::vector<int>& getAlphabetStates() const {return profileSet_[0]->getAlphabetStates(); }

    int getAlphabetStateAsInt(size_t i) const {return profileSet_[0]->getAlphabetStateAsInt(i); }

    std::string getAlphabetStateAsChar(size_t i) const {return profileSet_[0]->getAlphabetStateAsChar(i); }
    

    /* Parametrizable interface */
    
    void fireParameterChanged(const ParameterList& parameters);


  protected:


    /* Utility functions */
    
    void associateParameters_(SubstitutionModelSet * profile, size_t profileIndex);
  };

} //namespace bpp
  
#endif /*_NON_HOMOGENEOUS_SITE_PROFILE_MIXTURE_H_ */
