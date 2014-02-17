//============================================================================
// Name        : PlottingTools.cpp
// Author      : P Symonds
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "../interface/Objects/MET.h"
#include "../interface/Objects/Muon.h"
#include "../interface/Objects/Jets.h"
#include "../interface/Objects/CutFlow.h"
#include "../interface/Objects/Fit.h"

using namespace std;

int main() {
	//Have to make a dummy histogram so that everything works. Silly ROOT
	TH1F* dummy = new TH1F("dummy", "dummy", 10, 0, 1);
	delete dummy;

	AllSamples samples("central", "");

//	MET met;
//	met.allPlots(samples);

//	Jets jets;
//	jets.allPlots(samples);

//	Muon muon;
//	muon.setSelection("TTbar_plus_X_analysis/MuPlusJets/QCD non iso mu+jets ge3j/");
//	muon.allPlots(samples);

//	CutFlow cutflow;
//	cutflow.allPlots(samples);

//	do not do muon plots and fits at the same time as it messes up the eta distn.
	Fit fit;
	fit.allFits();

	return 0;

}
