/*
 * AllSamples.cpp
 *
 *  Created on: Oct 15, 2013
 *      Author: philip
 */

#include "../../interface/Samples/AllSamples.h"

namespace std {

AllSamples::AllSamples() {
	// TODO Auto-generated constructor stub

	single_mu_data = new Sample;
	Sample data_temp("SingleMu", kBlack, kBlack);
	*single_mu_data = data_temp;

	ttbar = new Sample;
	Sample ttbar_temp("TTJet", kRed, kBlack);
	*ttbar = ttbar_temp;

	single_t = new Sample;
	Sample single_t_temp("SingleTop", kMagenta, kBlack);
	*single_t = single_t_temp;

	vjets = new Sample;
	Sample vjets_temp("VJets", kGreen-3, kBlack);
	*vjets = vjets_temp;

	qcd = new Sample;
	Sample qcd_temp("QCD_Muon", kYellow, kBlack);
	*qcd = qcd_temp;
}

AllSamples::~AllSamples() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
