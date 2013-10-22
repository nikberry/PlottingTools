/*
 * MET.cpp
 *
 *  Created on: Oct 11, 2013
 *      Author: philip
 */

#include "../../interface/Objects/MET.h"
#include "TCanvas.h"
//#include "TLatex.h"
//#include <sstream>
//#include <iomanip>

namespace std {

MET::MET() {

}

MET::~MET() {

}

void MET::allPlots(AllSamples samples){


//  Variable::Variable(TString name_temp, TString xTitle_temp, double minX_temp, double maxX_temp, int rebinFact_temp)
	Variable met;
	savePlot(samples, met);

	Variable met_sig("METsignificance", "E_{T}^{miss} significance", 0, 150, 2);
	savePlot(samples, met_sig);
}

//void MET::savePlot(AllSamples samples, Variable variable){
//
//	TdrStyle style;
//	style.setTDRStyle();
//
//	//draw histos to files
//	TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);
//
//	TH1D* data = readGe2bHistogram(*samples.single_mu_data, variable);
//	samples.single_mu_data->SetHisto(data);
//
//	THStack *hs = buildStack(samples, variable);
//
//	hs->Draw();
//	data->Draw("E same");
//	data->SetMarkerStyle(20);
//
//	hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);
//	hs->GetXaxis()->SetLimits(variable.minX, variable.maxX);
//	hs->GetXaxis()->SetTitle(variable.xTitle); hs->GetXaxis()->SetTitleSize(0.05);
//	hs->GetYaxis()->SetTitle("Number of Events");hs->GetYaxis()->SetTitleSize(0.05);
//
//	TLegend* leg = legend(samples);
//	leg->Draw();
//
//	TText* textChan = doChan(0.12,0.96);
//	textChan->Draw();
//	TText* textPrelim = doPrelim(0.6,0.96);
//	textPrelim->Draw();
//
//	c1->SaveAs(variable.name+"_test.png");
//
//	delete hs;
//	delete c1;
//	delete leg;
//	delete textChan;
//	delete textPrelim;
//
//}


//TH1D* MET::readGe2bHistogram(Sample sample, Variable variable) {
//
////	cout << "plot: " << "TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable << endl;
//
//	TH1D* plot = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_2btags");
//	TH1D* plot2 = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_3btags");
//	TH1D* plot3 = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_4orMoreBtags");
//
//	plot->Add(plot2);
//	plot->Add(plot3);
//
//	plot->SetFillColor(sample.fillColor);
//	plot->SetLineColor(sample.lineColor);
//	plot->Rebin(variable.rebinFact);
//
//	return plot;
//}
//
//THStack* MET::buildStack(AllSamples samples, Variable variable){
//
//	THStack *hs = new THStack("hs","test");
//
//	TH1D* ttbar = readGe2bHistogram(*samples.ttbar, variable);
//	TH1D* single_t = readGe2bHistogram(*samples.single_t, variable);
//	TH1D* vjets = readGe2bHistogram(*samples.vjets, variable);
//	TH1D* qcd = readGe2bHistogram(*samples.qcd, variable);
//
//	samples.ttbar->SetHisto(ttbar);
//	samples.single_t->SetHisto(single_t);
//	samples.vjets->SetHisto(vjets);
//	samples.qcd->SetHisto(qcd);
//
//	hs->Add(qcd);
//	hs->Add(vjets);
//	hs->Add(single_t);
//	hs->Add(ttbar);
//
//	return hs;
//}



} /* namespace std */


