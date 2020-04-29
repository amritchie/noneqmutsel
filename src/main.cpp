#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include <getopt.h>
#include <unistd.h>

#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/GeneticCode.all>

#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Io.all>
#include <Bpp/Phyl/Distance.all>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Simulation.all>

#include <Bpp/Numeric/Prob.all>
#include <Bpp/App/ApplicationTools.h>

#include "/home/aritchie/src/nemss2/src/MutationSelectionModel.h"
#include "/home/aritchie/src/nemss2/src/AdditionalOptimizationTools.h"
//#include "/home/aritchie/src/nemss2/src/NonHomogeneousSiteProfileMixture.h"
#include "/home/aritchie/src/nemss2/src/RNonHomogeneousMultiProfileTreeLikelihood.h"
#include "/home/aritchie/src/nemss2/src/BranchMixtureMutSelInference.h"


#include "nemss_config.h"

using namespace std;
using namespace bpp;

/***********************************Procedure Division************************************/

void readBreakpointFile(string brkPath, vector <vector <int>* >& breakpointLists) {
  stringstream bps;
  ifstream brkFile;
  string brk;
  brkFile.open(brkPath);
  while (getline(brkFile, brk)) {
    vector<int>* breakpoints = new vector<int>();
    bps << brk;
    int bp;
    while (bps >> bp) {
      breakpoints->push_back(bp);
    }
    breakpointLists.push_back(breakpoints);
    bps.str("");
    bps.clear();
  }
  brkFile.close();
}

void splitSimInfoLine(string line, vector <string>& output) {

  string token;
  stringstream s(line);
  while (getline(s, token, ',')) {
    if (token.find('\r') != string::npos) token.pop_back(); 
    output.push_back(token);
  }

}

void readSimInfoFile(string path, vector<ParameterList * >& trueParams) {

  ifstream simInfoFile;
  simInfoFile.open(path);
    
  string simInfoLine;
    
  vector <string> headers;
  getline(simInfoFile, simInfoLine);
  splitSimInfoLine(simInfoLine, headers);

  while (getline(simInfoFile, simInfoLine)) {

    ParameterList * pl = new ParameterList();
    vector <int> * simCatSizes = new vector <int>;
    vector <string> simParams;
    splitSimInfoLine(simInfoLine, simParams);

    for (int i = 0; i < headers.size(); i++) {
      stringstream s;
      string pname = headers[i];
	
      auto fit = pname.find("fit");
      auto prob = pname.find("relproba");
      auto size = pname.find("size");
      
      if (fit != string::npos or prob != string::npos or size != string::npos) {
	double d = 0.0;
	s << simParams[i];
	s >> d;
	pl->addParameter(new Parameter(pname, d));
      }
	
    }

    trueParams.push_back(pl);
      
  }
    
  simInfoFile.close();
}

void readEstimatesFile(string path, ParameterList& estParams) {

  ifstream estimFile;
  estimFile.open(path);

  string line;
  getline(estimFile, line);
  getline(estimFile, line);

  while (getline(estimFile, line)) {
      string pname;
      double pval;
      
      stringstream s;
      s << line;
      s >> pname;
      s >> pval;

      estParams.addParameter(new Parameter(pname, pval));
    
  }

  estimFile.close();

}

void packParametersByModel(ParameterList allParams,
			   Tree * tree,
			   vector<ParameterList>& paramsByModel) {

  map < int, ParameterList > pList;
  
  for (int i=0; i < allParams.size(); i++) {
    
    Parameter p = allParams[i];
    string pname = p.getName();

    if (pname.find("fit") == string::npos)
      continue;

    int bpId = TextTools::toInt(pname.substr(2, pname.find("_") - 2));
    int submodId = TextTools:: toInt(pname.substr(pname.find("m") + 1,
						  pname.find("fit") - pname.find("m") - 2));
    string submodStr = pname.substr(pname.find("_") + 1, pname.find("fit") - pname.find("_") - 2); 

    if (bpId == tree->getRootId()) {

      if (pList.find(submodId) == pList.end()) {
	ParameterList newpl;
	pList[submodId] = newpl;
      }
      
      pList[submodId].addParameter(p);
      pList[submodId].getParameter(pname).setName(pname.replace(pname.find(submodStr),
								submodStr.size(),
								"1_m0"));
      
    }
    
  }

  for (auto it = pList.begin(); it != pList.end(); ++it) {
    paramsByModel.push_back(it->second);
  }

}

/*********************************** MAIN PROCEDURE **************************************/

/** Procedure finds likelihood with a breakpoint at each free branch, saves the position of 
 ** the most likely branch, then looks for the next likeliest branch to have a breakpoint.
 **
 ** Currently this just goes until we tell it to stop, then we can rank models with different 
 ** numbers of breakpoints via AICc. This is because we know how many BPs there are in simulated data;
 ** wanted is to build a penalised likelihood method to do this automatically.
 **/

int main (int argc, char* argv[]) {

  static string prefix;
  static string outputSuffix;

  static int treeIndexStart, treeIndexStop;

  static int independentRootFreqsFlag;
  static int inferBreakpointsFlag;
  static int inferFreqShiftsFlag;
  static int doNotRunFlag;

  static int numberOfCategories;
  static int maxBreakpoints;

  /* Parse arguments */

  /* ./src/mutsel --prefix test_1brk_2cat_unif --output validate --start 11 --end 11 --categories 2 --maxbps 1 */

  static struct option cmdOpts[] =
    {
      /* flags */
      {"rootfreqs", no_argument, &independentRootFreqsFlag, 1 },
      {"inferbps", optional_argument, 0, 'b'},
      {"inferfreqshifts", optional_argument, 0, 'f'},
      {"donotrun", no_argument, &doNotRunFlag, 1},
      /* arguments */
      {"prefix", required_argument, 0, 'p'},
      {"output", required_argument, 0, 'o'},
      {"start", required_argument, 0, 's'},
      {"end", required_argument, 0, 'e'},
      {"categories", required_argument, 0, 'k'},
      {"maxbps", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };

  int c;
  int optionIndex = 0;
  vector<int> assumedBpIndices;
  while (1) {
   
    c = getopt_long (argc, argv, "rb::dp:o:s:e:", cmdOpts, &optionIndex);
    if (c == -1) break;
    stringstream s;
    
    switch (c)
      {
      case 0:
	if (cmdOpts[optionIndex].flag != 0)
	  break;
      case 'p':
	prefix = optarg;
	break;
      case 'o':
	outputSuffix = optarg;
	outputSuffix = "_" + outputSuffix;
	break;
      case 's':
	s << optarg;
	s >> treeIndexStart;
	break;
      case 'e':
	s << optarg;
	s >> treeIndexStop;
	break;
      case 'b':
	if (optarg) {
	  s << optarg;
	  int tempBp;
	  while (s >> tempBp) {
	    assumedBpIndices.push_back(tempBp);
	  }
	}
	inferBreakpointsFlag=1;
	break;
      case 'f':
	if (optarg) {
	  s << optarg;
	  int tempBp;
	  while (s >> tempBp) {
	    assumedBpIndices.push_back(tempBp);
	  }
	}
	inferFreqShiftsFlag=1;
	break;
      case 'k':
	if(optarg) {
	  s << optarg;
	  s >> numberOfCategories;
	} else {
	  numberOfCategories = 1;
	}
	break;
      case 'm':
	if(optarg) {
	  s << optarg;
	  s >> maxBreakpoints;
	} else {
	  maxBreakpoints=0;
	}
	break;
      case '?':
	return -1;
      default:
	abort();
      }

    s.clear();
  }

  for (int index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);

  /* End argument parsing */


  try {

    string srcDir = string(get_current_dir_name()) + "/";
    string pathPrefix = srcDir + prefix;
    string treePath = pathPrefix + ".trees";

    /* Read trees and branch lengths from newick file */
    Newick treeReader;
    vector<Tree * > trees;
    treeReader.read(treePath, trees);

    /* Read breakpoints from associated _brk.csv file */
    vector<vector <int>* > breakpointLists;
    readBreakpointFile(pathPrefix + "_brk.csv", breakpointLists);

    /* Read true parameters from simulation info file */

    vector <ParameterList *> trueParamLists;
    readSimInfoFile(pathPrefix + "_siminfo.csv", trueParamLists);

    /* Running through trees and conducting inferences */
    for (int treeIndex=treeIndexStart; treeIndex<=treeIndexStop; treeIndex++) {

      Tree * tree = trees[treeIndex];
      vector <int> breakpoints = *(breakpointLists[treeIndex]);
      ParameterList trueParams = *(trueParamLists[treeIndex]);
      
      string treePrefix = pathPrefix + "_" + to_string(treeIndex);
      string alignPath = treePrefix + ".fasta";

      Fasta seqReader;
      CodonAlphabet  codonAlpha(&AlphabetTools::DNA_ALPHABET);
      AlignedSequenceContainer *sequences = seqReader.readAlignment(alignPath, &codonAlpha);
      SiteContainer *sites = new VectorSiteContainer(*sequences);
      delete sequences;

      // Update pyvolve parameter names to match BPP naming scheme

      int trueNBackgroundCats = 0;

      ParameterList renamedTrueParams;
      
      for (int p=0; p <trueParams.size(); p++) {
	
	string pname = trueParams[p].getName();
	string newParamName;
	auto isProb = pname.find("relproba");
	auto isFitness = pname.find("fit");

	if (isFitness != string::npos || isProb != string::npos) {

	  std::smatch sm;
	  std::regex bpRegex ("bp(\\d+)");

	  std::regex_search(pname, sm, bpRegex);
	  
	  int bpNumber = TextTools::toInt(sm[1].str());
	  int bpNodeNumber = bpNumber == 0?
	    tree->getNumberOfNodes() - 1 :
	    breakpoints[breakpoints.size() - bpNumber]; // Reverse postorder
	  
	  newParamName = std::regex_replace(pname, bpRegex, "bp" + TextTools::toString(bpNodeNumber));

	  if (isFitness != string::npos) {

	    vector<int> pyvolveToBpp( { 0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18 } );
	    
	    std::regex fitRegex ("fit(\\d+)");
	    std::regex_search(newParamName, sm, fitRegex);
	    
	    int fitState = TextTools::toInt(sm[1].str());
	    int bppFitState = pyvolveToBpp[fitState];

	    // Pyvolve state 0 is Alanine, which is fixed in BPP
	    if (bppFitState == 0)
	      continue;
	    
	    newParamName = std::regex_replace(newParamName,
					      fitRegex,
					      "fit" + TextTools::toString(bppFitState - 1));
	    
	  } else if (isProb != string::npos && bpNumber == 0) {
	 
	      trueNBackgroundCats++;
	      
	  }
	  
	  renamedTrueParams.addParameter( new Parameter(newParamName, trueParams[p].getValue()));
	  
	}

      }
      
      vector<ParameterList> paramsByModel;
      packParametersByModel(renamedTrueParams, tree, paramsByModel);

      vector<int> assumedBps;
      for (int i : assumedBpIndices) {
      	assumedBps.push_back(breakpoints[breakpoints.size() - i - 1]);
      }

      BranchMixtureMutSelInference inference(tree, sites, assumedBps);


      // Output individual site likelihoods
      ofstream siteLikelihoodsOut;

      try {
	siteLikelihoodsOut.open(treePrefix + "_siteLnLs" + outputSuffix + ".csv");
      } catch (Exception e) {
	cout << "Could not open analysis files with specified prefix." << endl;
      }


      VVdouble siteLnLs;
      siteLnLs.push_back(inference.computeSiteLogLikelihoodsForParameters(renamedTrueParams,
									  trueNBackgroundCats));
	
      for (ParameterList pl : paramsByModel) {
	siteLnLs.push_back(inference.computeSiteLogLikelihoodsForParameters(pl, 1));
      }

      siteLikelihoodsOut << "Full_model";
	
      for (int i=0; i < trueNBackgroundCats; i++) {
	siteLikelihoodsOut << ",Background_" << i;
      }
	
      siteLikelihoodsOut << endl;
      
      for(int i=0; i < sites->getNumberOfSites(); i++) {
	siteLikelihoodsOut << siteLnLs[0][i];
	for (int m=1; m < siteLnLs.size(); m++) {
	  siteLikelihoodsOut << "," << siteLnLs[m][i];
	}
	siteLikelihoodsOut << endl;
      }

      siteLikelihoodsOut.close();
	
      // If not running, we look at the empirical, true and estimated distributions at the tips
	
      if (doNotRunFlag) {

	 ofstream analyseOut;

	 try {
	   analyseOut.open(treePrefix + "_analyse" + outputSuffix + ".txt");
	 } catch (Exception e) {
	   cout << "Could not open analysis files with specified prefix." << endl;
	 }

	 double trueLnL = inference.computeLogLikelihoodForParameters(renamedTrueParams, trueNBackgroundCats);
	 analyseOut << "Log Likelihood of true parameters:" << endl;
	 analyseOut << trueLnL << endl << endl;
	
	 ParameterList prevEstParams;
	 VVdouble estTipDistrs;
	 readEstimatesFile(treePrefix + "_estim" + outputSuffix + ".txt", prevEstParams);
	 inference.computeTipDistrsForParameters(prevEstParams, numberOfCategories, estTipDistrs);
	 VVdouble trueTipDistrs;
	 inference.computeTipDistrsForParameters(renamedTrueParams,
						 trueNBackgroundCats,
						 trueTipDistrs);

	 VVdouble empTipDistrs;
	 Vdouble avTipDistrs;
	 inference.getEmpiricalCodonFreqs(empTipDistrs, avTipDistrs);
	
	 analyseOut << "AA distributions for individual tips (based on previous estimates):" << endl;
	
	 for (int i=0; i < estTipDistrs.size(); i++) {
	   analyseOut << i;
	   double kld = 0;
	   double sum = 0;
	   for (int j=0; j < estTipDistrs[i].size(); j++) {
	     analyseOut << "\t" << estTipDistrs[i][j];
	     if (estTipDistrs[i][j] > NumConstants::VERY_TINY()
	 	&& trueTipDistrs[i][j] > NumConstants::VERY_TINY()) {
	       kld += -estTipDistrs[i][j] * (log(trueTipDistrs[i][j]) - log(estTipDistrs[i][j]));
	     }
	     sum += trueTipDistrs[i][j];
	   }
	   analyseOut << "\t" << kld << endl;
	 }

	 analyseOut << "AA distributions for individual tips (based on true parameters):" << endl;

	 for (int i=0; i < trueTipDistrs.size(); i++) {
	   analyseOut << i;
	   double kld = 0;
	   for (int j=0; j < trueTipDistrs[i].size(); j++) {
	     analyseOut << "\t" << trueTipDistrs[i][j];
	     if (estTipDistrs[i][j] > NumConstants::VERY_TINY()
	 	&& trueTipDistrs[i][j] > NumConstants::VERY_TINY()) {
	       kld += -trueTipDistrs[i][j] * (log(estTipDistrs[i][j]/trueTipDistrs[i][j]));
	     }
	   }
	   analyseOut << "\t" << kld << endl;
	 }

	 analyseOut << "Average empirical codon frequencies across all tips: " << endl;
	
	 for (int i=0; i < avTipDistrs.size(); i++) {
	   analyseOut << "\t" << avTipDistrs[i];
	 }

	 analyseOut.close();
	
      } else { // we are running inferences

	if (inferBreakpointsFlag) { // infer breakpoints and parameters, currently no mixtures.
	
	  inference.setAssumedBreakpoints(assumedBps);
	  inference.bppTestBestBranch(maxBreakpoints, 1, 1);
	  inference.writeBreakpointEstimates(treePrefix + "_estbps" + outputSuffix + ".txt");
	
	} else if (inferFreqShiftsFlag) { // infer mixture frequency shifts with known mixture

	  inference.setAssumedBreakpoints(assumedBps);
	  inference.bppFindSiteFrequencyShifts(renamedTrueParams,
					       trueNBackgroundCats,
					       maxBreakpoints);
	  inference.writeBreakpointEstimates(treePrefix + "_estbps" + outputSuffix + ".txt");

	} else { // If nothing else flagged, infer the background with no breakpoints

	  inference.setAssumedBreakpoints(assumedBps);
	  inference.bppInferBackgroundMixture(numberOfCategories, independentRootFreqsFlag);
	  
	}

	inference.writeParameterEstimates(treePrefix + "_estim" + outputSuffix + ".txt");
	inference.writeEquilibriumFrequencies(treePrefix + "_estfreqs" + outputSuffix + ".txt");
	inference.writeTimings(treePrefix + "_timings" + outputSuffix + ".csv");
      
	}
    }
    
  }
  catch(Exception& e) {
    cout << "Bio++ exception:" << endl;
    cout << e.what() << endl;
    return(-1);
  }
  catch(exception e) {
    cout << "Some other exception:" << endl;
    cout << e.what() << endl;
    return(-1);
  }
  
  return 0;
 
}

/************************************* END OF FILE ***************************************/
