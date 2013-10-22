/*
 * getSample.cpp
 *
 *  Created on: Oct 14, 2013
 *      Author: philip
 */

#include "../../interface/Samples/Sample.h"


namespace std {

Sample::Sample() {

	file = new TFile;// ("~/workspace/TtbarWork/rootFilesV4/central/TTJet_19584pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	histo = new TH1D;
	fillColor = kRed;
	lineColor = kBlack;
}

Sample::Sample(TString sample, Color_t fColor, Color_t lColor) {
	cout << "Getting file: " << sample + "_19584pb_PFElectron_PFMuon_PF2PATJets_PFMET.root" << endl;
	file = new TFile ("~/workspace/TtbarWork/rootFilesV4/central/" + sample + "_19584pb_PFElectron_PFMuon_PF2PATJets_PFMET.root");
	histo = new TH1D;
	fillColor = fColor;
	lineColor = lColor;
}

void Sample::SetHisto(TH1D* plot){
	histo = plot;
}

void Sample::SetFillColor(Color_t fColor){
	fillColor = fColor;
}

Color_t Sample::GetFillColor(){
	return fillColor;
}

void Sample::SetLineColor(Color_t lColor){
	lineColor = lColor;
}

Color_t Sample::GetLineColor(){
	return lineColor;
}

Sample::~Sample() {

}


} /* namespace std */
