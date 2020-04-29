/** File: NonHomogeneousSiteProfileMixture **/

#include <string>
#include <vector>
#include <set>

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Text/TextTools.h>

#include <Bpp/Phyl/Model/SubstitutionModelSet.h>

#include "/home/aritchie/src/nemss2/src/NonHomogeneousSiteProfileMixture.h"


using namespace bpp;

using namespace std;


NonHomogeneousSiteProfileMixture::NonHomogeneousSiteProfileMixture(vector <SubstitutionModelSet*>& profileSet) :
  AbstractParameterAliasable(""),
  profileSet_(),
  profileParameters_()
{

  string alphaType = profileSet[0]->getAlphabet()->getAlphabetType();
  
  for (auto p : profileSet) {
    if (p->getAlphabet()->getAlphabetType() != alphaType)
      throw
	Exception ("NonHomogeneousSiteProfileMixture::NonHomogeneousSiteProfileMixture.All profiles in mixture must use the same alphabet type.");
  }
  
  profileParameters_.resize(profileSet.size());

  for (size_t i=0; i<profileSet.size(); i++) addProfile(profileSet[i]);
}


NonHomogeneousSiteProfileMixture::NonHomogeneousSiteProfileMixture(const NonHomogeneousSiteProfileMixture& mix) :
AbstractParameterAliasable(mix),
  profileSet_(),
  profileParameters_()
{

  profileParameters_.resize(mix.profileSet_.size());
  
  for (size_t i=0; i< mix.getNumberOfProfiles(); i++) {
    profileSet_.push_back(mix.profileSet_[i]->clone());
    profileParameters_[i] = *(mix.profileSet_[i]->getIndependentParameters().clone());
  }

}


void NonHomogeneousSiteProfileMixture::fireParameterChanged(const ParameterList& parameters) {

  AbstractParameterAliasable::fireParameterChanged(parameters);

  set <int> catsToUpdate;

  vector <string> names = parameters.getParameterNames();
  names = VectorTools::vectorUnion<string>(names, getAliasedParameters(parameters).getParameterNames());
  
  for (string n : names) {

    string pName = n.substr(3, n.npos);
    int cat = stoi(n.substr(1));
    profileParameters_[cat].getParameter(pName).setValue(getParameterValue(n));
    catsToUpdate.insert(cat);
    
  }

  for (int k : catsToUpdate) {
    
    profileSet_[k]->matchParametersValues(profileParameters_[k]);
    
  }
  
}


void NonHomogeneousSiteProfileMixture::addProfile(SubstitutionModelSet * profile) {

  profileSet_.push_back(profile);

  associateParameters_(profile, profileSet_.size() - 1);

}


void NonHomogeneousSiteProfileMixture::associateParameters_(SubstitutionModelSet * profile,
							    size_t profileIndex) {

  ParameterList pl = profile->getIndependentParameters();

  profileParameters_[profileIndex] = *(pl.clone());
  
  vector <string> names = pl.getParameterNames();
  
  for (auto n : names) {
   
    Parameter * newP = new Parameter(pl.getParameter(n));
    newP->setName("p" + TextTools::toString(profileIndex) + "_" + n);
    addParameter_(newP);
    
  }

}


vector <int>  NonHomogeneousSiteProfileMixture::getNodesWithParameter(string pname) const {

  if (!hasIndependentParameter(pname))
    throw
      Exception("NonHomogeneousSiteProfileMixture::getNodesWithParameter(string pname). No parameter with given name.");

  vector <string> aliases;
  aliases.push_back(pname);
  aliases = VectorTools::vectorUnion<string>(aliases, getAlias(pname));

  vector<int> nodeIds;
  for (string al : aliases) {

    string subName = al.substr(3, pname.npos);
    size_t cat = stoi(al.substr(1));
    nodeIds = VectorTools::vectorUnion<int>(nodeIds, profileSet_[cat]->getNodesWithParameter(subName));
    
  }

  return nodeIds;
}


const TransitionModel * NonHomogeneousSiteProfileMixture::getModelForProfileForNode(size_t profile, int nodeId) const {

  return profileSet_[profile]->getModelForNode(nodeId);

}


const SubstitutionModel * NonHomogeneousSiteProfileMixture::getSubstitutionModelForProfileForNode(size_t profile, int nodeId) const {

  return profileSet_[profile]->getSubstitutionModelForNode(nodeId);

}


TransitionModel * NonHomogeneousSiteProfileMixture::getModelForProfileForNode(size_t profile, int nodeId) {

  return profileSet_[profile]->getModelForNode(nodeId);

}


SubstitutionModel * NonHomogeneousSiteProfileMixture::getSubstitutionModelForProfileForNode(size_t profile, int nodeId) {

  return profileSet_[profile]->getSubstitutionModelForNode(nodeId);

}


vector <TransitionModel * > NonHomogeneousSiteProfileMixture::getAllModelsForNode(int nodeId) {

  vector <TransitionModel * > modelsForNode;
  
  for (auto p : profileSet_) {
    modelsForNode.push_back(p->getModelForNode(nodeId));
  }

  return modelsForNode;

}


vector <SubstitutionModel * > NonHomogeneousSiteProfileMixture::getAllSubstitutionModelsForNode(int nodeId){

  vector <SubstitutionModel * > subModelsForNode;

  for (auto p : profileSet_) {
    subModelsForNode.push_back(p->getSubstitutionModelForNode(nodeId));
  }

  return subModelsForNode;
  
}


ParameterList NonHomogeneousSiteProfileMixture::getNodeParametersForProfile (size_t profile) const {

  ParameterList pNodeParams = profileSet_[profile]->getNodeParameters();
  ParameterList aliasedFrom = getFromParameters(pNodeParams);
 
  ParameterList pl;
  
  vector <string> names = VectorTools::vectorUnion<string>(pNodeParams.getParameterNames(),
							   aliasedFrom.getParameterNames());
  
  for (auto n : names) {
   
    string pname = "p" + TextTools::toString(profile) + "_" + n;

    if (hasIndependentParameter(pname)) {
      Parameter * newP = new Parameter(pNodeParams.getParameter(n));
      newP->setName(pname);
      pl.addParameter(newP);
    }
    
  }

  return pl;

}


const vector <double> NonHomogeneousSiteProfileMixture::getRootFrequenciesForProfile (size_t profile) const {

  return  profileSet_[profile]->getRootFrequencies();

}

/** 
 * Convenience function - where all profiles have a parameter with the same underlying name,
 * tie them all together so they act as one global parameter.
 *
 */
void NonHomogeneousSiteProfileMixture::aliasParameterAcrossProfiles(string name) {

  string globalParName = "p0_" + name;
  for (int k=1; k < getNumberOfProfiles(); k++) {
    string localParName = "p" + TextTools::toString(k) + "_" + name;
    if (hasParameter(localParName)) {
	aliasParameters(globalParName, localParName);
      } else {
      throw Exception("NonHomogeneousSiteProfileMixture::aliasParameterAcrossProfiles. Parameter name must be shared across all profiles in mixture.");
    }
  }

}

bool NonHomogeneousSiteProfileMixture::allAreFullySetUpFor(const Tree& tree) const {

  for (auto p : profileSet_) {
    if (!p->isFullySetUpFor(tree)) return false;
  }

  return true;
  
}
