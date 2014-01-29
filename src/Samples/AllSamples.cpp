/*
 * AllSamples.cpp
 *
 *  Created on: Oct 15, 2013
 *      Author: philip
 */

#include "../../interface/Samples/AllSamples.h"

namespace std {

AllSamples::AllSamples(TString systematic, TString eSystematic) {

	single_mu_data = new Sample;
	Sample data_temp("SingleMu", kBlack, kBlack, "central", "");
	*single_mu_data = data_temp;

	ttbar = new Sample;
	Sample ttbar_temp("TTJet", kRed, kBlack, systematic, eSystematic);
	*ttbar = ttbar_temp;

	single_t = new Sample;
	Sample single_t_temp("SingleTop", kMagenta, kBlack, systematic, eSystematic);
	*single_t = single_t_temp;

	signal = new Sample;
	Sample signal_temp("TTJet", kRed, kBlack, systematic, eSystematic);
	*ttbar = signal_temp;

	vjets = new Sample;
	Sample vjets_temp("VJets", kGreen-3, kBlack, systematic, eSystematic);
	*vjets = vjets_temp;

	qcd = new Sample;
	Sample qcd_temp("QCD_Muon", kYellow, kBlack, systematic, eSystematic);
	*qcd = qcd_temp;
}

AllSamples::~AllSamples() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
