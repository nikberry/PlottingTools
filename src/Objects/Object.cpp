/*
 * Object.cpp
 *
 *  Created on: Oct 18, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Object.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <sstream>
#include <iomanip>

namespace std {

Object::Object() {
	// TODO Auto-generated constructor stub

}

Object::~Object() {
	// TODO Auto-generated destructor stub
}

TH1D* Object::readGe2bHistogram(Sample sample, Variable variable) {

	TH1D* plot = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_2btags");
	TH1D* plot2 = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_3btags");
	TH1D* plot3 = (TH1D*) sample.file->Get("TTbar_plus_X_analysis/MuPlusJets/Ref selection/MET/patType1CorrectedPFMet/"+variable.name+"_4orMoreBtags");

	plot->Add(plot2);
	plot->Add(plot3);

	plot->SetFillColor(sample.fillColor);
	plot->SetLineColor(sample.lineColor);
	plot->Rebin(variable.rebinFact);

	return plot;
}

THStack* Object::buildStack(AllSamples samples, Variable variable){

	THStack *hs = new THStack("hs","test");

	TH1D* ttbar = readGe2bHistogram(*samples.ttbar, variable);
	TH1D* single_t = readGe2bHistogram(*samples.single_t, variable);
	TH1D* vjets = readGe2bHistogram(*samples.vjets, variable);
	TH1D* qcd = readGe2bHistogram(*samples.qcd, variable);

	samples.ttbar->SetHisto(ttbar);
	samples.single_t->SetHisto(single_t);
	samples.vjets->SetHisto(vjets);
	samples.qcd->SetHisto(qcd);

	hs->Add(qcd);
	hs->Add(vjets);
	hs->Add(single_t);
	hs->Add(ttbar);

	return hs;
}

void Object::savePlot(AllSamples samples, Variable variable){

	TdrStyle style;
	style.setTDRStyle();

	//draw histos to files
	TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

	TH1D* data = readGe2bHistogram(*samples.single_mu_data, variable);
	samples.single_mu_data->SetHisto(data);

	THStack *hs = buildStack(samples, variable);

	hs->Draw();
	data->Draw("E same");
	data->SetMarkerStyle(20);

	hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);
	hs->GetXaxis()->SetLimits(variable.minX, variable.maxX);
	hs->GetXaxis()->SetTitle(variable.xTitle); hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitle("Number of Events");hs->GetYaxis()->SetTitleSize(0.05);

	TLegend* leg = legend(samples);
	leg->Draw();

	TText* textChan = doChan(0.12,0.96);
	textChan->Draw();
	TText* textPrelim = doPrelim(0.6,0.96);
	textPrelim->Draw();

	c1->SaveAs(variable.name+"_test.png");

	delete hs;
	delete c1;
	delete leg;
	delete textChan;
	delete textPrelim;

}

TLegend* Object::legend(AllSamples samples){

		TLegend *tleg;
		tleg = new TLegend(0.7,0.7,0.8,0.9);
		tleg->SetTextSize(0.04);
		tleg->SetBorderSize(0);
		tleg->SetFillColor(10);
		tleg->AddEntry(samples.single_mu_data->histo , "2012 data", "lpe");
		tleg->AddEntry(samples.ttbar->histo , "t#bar{t}", "f");
		tleg->AddEntry(samples.single_t->histo, "single top"      , "f");
		tleg->AddEntry(samples.vjets->histo , "v+jets", "f");
		tleg->AddEntry(samples.qcd->histo, "QCD"      , "f");

		return tleg;
}

TText* Object::doChan(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "#mu, #geq 4 jets, #geq 2 btags";

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

TText* Object::doPrelim(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "CMS Preliminary, L = 19.6 fb^{-1}";

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

} /* namespace std */
