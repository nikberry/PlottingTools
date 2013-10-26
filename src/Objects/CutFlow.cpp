/*
 * CutFlow.cpp
 *
 *  Created on: Oct 26, 2013
 *      Author: philip
 */

#include "../../interface/Objects/CutFlow.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include <sstream>
#include <iomanip>
#include <math.h>

namespace std {

CutFlow::CutFlow() {
	// TODO Auto-generated constructor stub
	objName = "CutFlow";
	selection = "EventCount";

}

CutFlow::~CutFlow() {
	// TODO Auto-generated destructor stub
}

void CutFlow::allPlots(AllSamples samples){

//  Variable::Variable(TString name_temp, TString xTitle_temp, double minX_temp, double maxX_temp, int rebinFact_temp)
	Variable ttbar("TTbarMuPlusJetsRefSelection", "selection step", -1.5, 9.5, 1);
	saveCutFlowPlot(samples, ttbar);

}

void CutFlow::saveCutFlowPlot(AllSamples samples, Variable variable){

	readCutFlowHistos(samples, variable);

	TH1D* data = samples.single_mu_data->histo;
	THStack *hs = buildStack(samples, variable);

	standardCutFlowPlot(data, hs, samples, variable);

	if(Globals::addRatioPlot){
		ratioCutFlowPlot(data, hs, samples, variable);
	}

	delete data;
	delete hs;
}

TH1D* CutFlow::readCutFlowHistogram(Sample sample, Variable variable) {

	cout << "plot: " << selection+"/"+variable.name << endl;

	TH1D* plot = (TH1D*) sample.file->Get(selection+"/"+variable.name);

	plot->SetFillColor(sample.fillColor);
	plot->SetLineColor(sample.lineColor);

	return plot;
}

void CutFlow::readCutFlowHistos(AllSamples samples, Variable variable){

	TH1D* data = readCutFlowHistogram(*samples.single_mu_data, variable);
	TH1D* ttbar = readCutFlowHistogram(*samples.ttbar, variable);
	TH1D* single_t = readCutFlowHistogram(*samples.single_t, variable);
	TH1D* vjets = readCutFlowHistogram(*samples.vjets, variable);
	TH1D* qcd = readCutFlowHistogram(*samples.qcd, variable);

	samples.single_mu_data->SetHisto(data);
	samples.ttbar->SetHisto(ttbar);
	samples.single_t->SetHisto(single_t);
	samples.vjets->SetHisto(vjets);
	samples.qcd->SetHisto(qcd);
}

THStack* CutFlow::buildStack(AllSamples samples, Variable variable){

	THStack *hs = new THStack("hs","test");

	hs->Add(samples.qcd->histo);
	hs->Add(samples.vjets->histo);
	hs->Add(samples.single_t->histo);
	hs->Add(samples.ttbar->histo);

	return hs;
}

TH1D* CutFlow::allMChisto(AllSamples samples, Variable variable){

	TH1D *allMC = (TH1D*)samples.ttbar->histo->Clone("ratio plot");

	allMC->Add(samples.qcd->histo);
	allMC->Add(samples.vjets->histo);
	allMC->Add(samples.single_t->histo);

	return allMC;
}

void CutFlow::standardCutFlowPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable){
	//Style
	TdrStyle style;
	style.setTDRStyle();

	//draw histos to files
	TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

	data->Draw();
	hs->Draw("hist");

	setBinLabels(hs, data);

	if(Globals::addHashErrors){
		TH1D* hashErrs = hashErrors(samples, variable);
		hashErrs->Draw("same e2");
	}

	data->Draw("E same");
	data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);

	hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);
	hs->GetXaxis()->SetLimits(variable.minX, variable.maxX);
	hs->GetXaxis()->SetTitle(variable.xTitle); hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitle("Number of Events");hs->GetYaxis()->SetTitleSize(0.05);

	TLegend* leg = legend(samples);
	leg->Draw();

	TText* textChan = doChan(0.12,0.96);
	textChan->Draw();
	TText* textPrelim = doPrelim(0.58,0.96);
	textPrelim->Draw();

	if(Globals::doLogPlot){
		c1->SetLogy();
		c1->SaveAs("Plots/ControlPlots/"+objName+"/Log/"+variable.name+".png");
		c1->SaveAs("Plots/ControlPlots/"+objName+"/Log/"+variable.name+".pdf");
	}else{
		c1->SaveAs("Plots/ControlPlots/"+objName+"/"+variable.name+".png");
		c1->SaveAs("Plots/ControlPlots/"+objName+"/"+variable.name+".pdf");
	}

	delete c1;
	delete leg;
	delete textChan;
	delete textPrelim;
}

void CutFlow::ratioCutFlowPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable){
	//draw histos with ratio plot
	float r = 0.3;
	float epsilon = 0.02;
	TCanvas *c2 = new TCanvas("Plot","Plot",635, 600);
	c2->SetFillColor(0);
	c2->SetFrameFillStyle(0);
	TPad *pad1 = new TPad("pad1","pad1",0,r-epsilon,1,1);
	pad1->SetBottomMargin(epsilon);
	c2->cd();
	pad1->Draw();
	pad1->cd();

	hs->Draw("hist");

	hs->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.3);
	hs->GetXaxis()->SetLimits(variable.minX, variable.maxX);
	hs->GetXaxis()->SetTitle(variable.xTitle); hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitle("Number of Events");hs->GetYaxis()->SetTitleSize(0.05);

	if(Globals::addHashErrors){
		TH1D* hashErrs = hashErrors(samples, variable);
		hashErrs->Draw("same e2");
	}

	data->Draw("E same");
	//data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);

	TLegend* leg = legend(samples);
	leg->Draw();

	TText* textChan = doChan(0.12,0.96);
	textChan->Draw();
	TText* textPrelim = doPrelim(0.58,0.96);
	textPrelim->Draw();

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,r*(1-epsilon));
	pad2->SetTopMargin(0);
	pad2->SetFrameFillStyle(4000);
	pad2->SetBottomMargin(0.4);
	c2->cd();
	pad2->Draw();
	pad2->cd();

	TH1D * allMC = allMChisto(samples, variable);
	TH1D * ratio = (TH1D*)data->Clone("ratio plot");
	ratio->Sumw2();
	ratio->SetStats(0);
	ratio->Divide(allMC);

	ratio->SetMaximum(2);
	ratio->SetMinimum(0.);

	setBinLabels(hs, ratio);

//	Will need to see if this works in other situations

	ratio->SetLabelSize(0.1, "X");
	ratio->SetTitleOffset(0.5, "Y");
	ratio->SetTitleOffset(0.8, "X");
	ratio->GetYaxis()->SetTitle("data/MC");ratio->GetYaxis()->SetTitleSize(0.1);
	ratio->GetXaxis()->SetTitle(variable.xTitle);ratio->GetXaxis()->SetTitleSize(0.15);

	ratio->Draw("ep");

	TLine *line = new TLine(variable.minX,1,variable.maxX,1);
	line->Draw();

	pad1->cd();

	if(Globals::doLogPlot){
		pad1->SetLogy();
		c2->SaveAs("Plots/ControlPlots/"+objName+"/Log/"+variable.name+"_ratio.png");
		c2->SaveAs("Plots/ControlPlots/"+objName+"/Log/"+variable.name+"_ratio.pdf");
	}else{
		c2->SaveAs("Plots/ControlPlots/"+objName+"/"+variable.name+"_ratio.png");
		c2->SaveAs("Plots/ControlPlots/"+objName+"/"+variable.name+"_ratio.pdf");
	}

	delete c2;
	delete leg;
	delete textChan;
	delete textPrelim;
}

TH1D* CutFlow::hashErrors(AllSamples samples, Variable variable){
	TH1D * hashErrors = allMChisto(samples, variable);

	hashErrors->SetFillColor(kBlack);
	hashErrors->SetFillStyle(3354);
	hashErrors->SetMarkerSize(0.);
	hashErrors->SetStats(0);

	return hashErrors;
}

TLegend* CutFlow::legend(AllSamples samples){

		TLegend *tleg;
		tleg = new TLegend(0.75,0.75,0.85,0.9);
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

TText* CutFlow::doChan(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "#mu, #geq 4 jets, #geq 2 btags";

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

TText* CutFlow::doPrelim(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "CMS Preliminary, L = "+Globals::lumi;

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

void CutFlow::setSelection(TString sel_name){
	selection = sel_name;
}

void CutFlow::setBinLabels(THStack* hs, TH1D* data){
	TString step[11] = {"Skim" ,"Cleaning and HLT","one isolated #mu", "loose #mu veto", "loose e veto", "#geq 1 jets", "#geq 2 jets","#geq 3 jets", "#geq 4 jets", "#geq1 CSV b-tag", "#geq2 CSV b-tag" };
	for(int i =0; i<data->GetNbinsX(); i++){
		hs->GetXaxis()->SetBinLabel(i+1, step[i]);
		data->GetXaxis()->SetBinLabel(i+1, step[i]);
	}

}

} /* namespace std */
