#ifndef SAMPLE_STEP2_H
#define SAMPLE_STEP2_H
#include "TFile.h"
#include "TChain.h"
#include "TTH/MEAnalysis/interface/METree.hh"
#include "TTH/Plotting/interface/easylogging++.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTH/TTHNtupleAnalyzer/interface/tth_tree.hh"

enum SampleType {
    NOME_8TEV,
    ME_8TEV,
    NOME_13TEV,
    ME_13TEV
};

enum Process {
    TTHBB,
    TTJETS
};

using namespace std;

class Sample {
public:
	const string nickName;
	const vector<string> fileNamesS1;
	const vector<string> fileNamesS2;
	const double fractionToProcess;
	long long totalEvents;
	const SampleType type;
	const Process process;
	const bool skip;

	const bool step1Enabled;
	const bool step2Enabled;

	TChain *chainS1, *chainS2;

	METree* treeS2;

	Sample(const edm::ParameterSet& pars);

	long long getFirstEvent() {
		if(skip) {
			return 0;
		}
		return 0;
	}
	long long getLastEvent() {
		if(skip) {
			return -1;
		}
		if (fractionToProcess < 1.0) {
			LOG(INFO) << "Processing a fraction of sample " << nickName << ": " << fractionToProcess; 
			return (long long)(fractionToProcess * totalEvents);
		}
		return totalEvents;
	}
};


#endif