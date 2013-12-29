/*
 * Fit.cpp
 *
 *  Created on: Dec 24, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Fit.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TMap.h"
#include "TVirtualFitter.h"
#include "TFitterMinuit.h"
#include <sstream>
#include <iomanip>
#include <math.h>

namespace std {


void fcn(int& npar, double* deriv, double& f, double par[], int flag){

		TObjArray* fithistos =  static_cast<TObjArray*> (TVirtualFitter::GetFitter()->GetObjectFit());

		double lnL = 0.0;

		TH1D* data = static_cast<TH1D*>(fithistos->At(0));
	    TH1D* top_fit = static_cast<TH1D*>(fithistos->At(1));
	    TH1D* bg_fit = static_cast<TH1D*>(fithistos->At(2));
	    TH1D* qcd_fit = static_cast<TH1D*>(fithistos->At(3));

	    double Nsignal = top_fit->Integral();
	    double Nvjets = bg_fit->Integral();
	    double Nqcd = qcd_fit->Integral();

	  for (int i=0; i< data->GetNbinsX(); i++){
	    //data_i is the observed number of events in each bin
	    int data_i = data->GetBinContent(i+1);
	    //xi is the expected number of events in each bin

	    double xi = par[0]*(top_fit->GetBinContent(i+1)/Nsignal) + par[1]*(bg_fit->GetBinContent(i+1)/Nvjets)
	    		+par[2]*(qcd_fit->GetBinContent(i+1)/Nqcd);


	    if(data_i !=0 && xi != 0){
	      lnL += log(TMath::Poisson(data_i, xi));
	    }

	  }

	  f = -2.0 * lnL;

	  //constraints
	   double nvjets_err =Nvjets*0.5;
	   double nqcd_err = Nqcd*2.;

//	cout << "nsig: " << Nsignal << " ,par0: " << par[0] << endl;
//	cout << "nvjets: " << Nvjets << " ,par1: " << par[1] << endl;
//	cout << "nqcd: " << Nqcd << " ,par2: " << par[2] << endl;

	  f += ((par[2]-Nqcd)*(par[2]-Nqcd))/(nqcd_err*nqcd_err);
	  f += ((par[1]-Nvjets)*(par[1]-Nvjets))/(nvjets_err*nvjets_err);

}

Fit::Fit() {
	objName = "Muon";
	selection = "TTbar_plus_X_analysis/MuPlusJets/Ref selection/";
	folder = "Fits";
}

Fit::~Fit() {

}


void Fit::allPlots(AllSamples samples){
	Variable absEta("muon_AbsEta", "muon |#eta|", 0, 2.6, 10);
	//standard plot without fit
	savePlot(samples, absEta);

}

void Fit::savePlot(AllSamples samples, Variable variable){

	readHistos(samples, variable);

	TH1D* data = samples.single_mu_data->histo;
	THStack *hs = buildStack(samples, variable);

	standardPlot(data, hs, samples, variable);

	if(Globals::addRatioPlot){
		ratioPlot(data, hs, samples, variable);
	}

	doFit(samples, variable);

	delete data;
	delete hs;
}

void Fit::doFit(AllSamples samples, Variable variable){

	//readHistos(samples, variable);
	//draw the templates used in the fit
	drawTemplates(samples, variable);

	TH1D* data = samples.single_mu_data->histo;
	TH1D* vjets = samples.vjets->histo;
	TH1D* qcd =   samples.qcd->histo;
	TH1D* signal = samples.signal->histo;

	//set the parameters
	double Ntop_err, Nbg_err, Nqcd_err;
	double Ntop = signal->IntegralAndError(0, signal->GetNbinsX()+1, Ntop_err);
	double Nbg  = vjets->IntegralAndError(0, signal->GetNbinsX()+1, Nbg_err);
	double Nqcd = qcd->IntegralAndError(0, signal->GetNbinsX()+1, Nqcd_err);
	double Ntotal = data->Integral();

	TVirtualFitter* fitter = TVirtualFitter::Fitter(0,3);
	fitter->SetFCN(fcn);

//	Double_t arg(-1); // disable printout
//	fitter->ExecuteCommand("SET PRINT",&arg,1);

	fitter->SetParameter(0,"Ntop", Ntop, Ntop_err,0,Ntotal);
	fitter->SetParameter(1,"Nbg" , Nbg,  Nbg_err,0,Ntotal);
	fitter->SetParameter(2,"Nqcd", Nqcd, Nqcd_err,0,Ntotal);

	TH1D* fit_histos[4] = {data, signal, vjets, qcd};
	TObjArray fithists(0);
	for(int i =0; i<4; i++){
		fithists.Add(fit_histos[i]);
	}

	//so that can use stuff in fcn
	fitter->SetObjectFit(&fithists);

	Double_t arglist[10];
	arglist[0] = 2;

//	fitter->FixParameter(2);
//	fitter->ExecuteCommand("CALL FCN", arglist, 1);
	fitter->ExecuteCommand("MIGRAD", arglist,0);
//	fitter->ExecuteCommand("MINOS", arglist, 0);

	 double results[] = { fitter->GetParameter(0), fitter->GetParameter(1),  fitter->GetParameter(2)};
	 double errors[] = { fitter->GetParError(0), fitter->GetParError(1), fitter->GetParError(2)};

	  double amin, edm, errdef;
	  int nvpar, nparx;

	  fitter->GetStats(amin, edm, errdef, nvpar, nparx);

	  //do the plotting bit
	  	signal->Scale(results[0]/Ntop);
		vjets->Scale(results[1]/Nbg);
		qcd->Scale(results[2]/Nqcd);

		TCanvas *c3 = new TCanvas("Plot","Plot",900, 600);

		  THStack* sum_fit = new THStack("sum fit","stacked histograms"); //used for stack plot
		  sum_fit->Add(qcd);sum_fit->Add(vjets);sum_fit->Add(signal);

		  sum_fit->Draw("hist");
		  data->Draw("E same");

		  sum_fit->GetXaxis()->SetLimits(variable.minX, variable.maxX);
		  sum_fit->GetXaxis()->SetTitle(variable.xTitle); sum_fit->GetXaxis()->SetTitleSize(0.05);
		  sum_fit->GetYaxis()->SetTitle("Number of Events");sum_fit->GetYaxis()->SetTitleSize(0.05);

			TLegend* leg = legend(samples);
			leg->Draw();

		 	TText* textChan = doChan(0.12,0.96);
			textChan->Draw();
			TText* textPrelim = doPrelim(0.58,0.96);
			textPrelim->Draw();

		 	c3->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_fit.pdf");
		    c3->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_fit.png");
		    delete c3;



}


void Fit::drawTemplates(AllSamples samples, Variable variable){
	TH1D* signal = (TH1D*)samples.signal->histo->Clone("signal");
	TH1D* vjets = (TH1D*)samples.vjets->histo->Clone("vjets");
	TH1D* qcd = (TH1D*)samples.qcd->histo->Clone("qcd");

	normAndColor(signal, *samples.ttbar);
	normAndColor(vjets, *samples.vjets);
	normAndColor(qcd, *samples.qcd);

	TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

	  signal->Draw();
	  vjets->Draw("same");
	  qcd->Draw("same");

	  signal->SetAxisRange(variable.minX, variable.maxX);
	  signal->GetXaxis()->SetTitle(variable.xTitle); signal->GetXaxis()->SetTitleSize(0.05);
	  signal->GetYaxis()->SetTitle("Normalised Events");signal->GetYaxis()->SetTitleSize(0.05);

	  	TLegend *tleg;
		tleg = new TLegend(0.65,0.7,0.8,0.9);
		tleg->SetTextSize(0.04);
		tleg->SetBorderSize(0);
		tleg->SetFillColor(10);

		tleg->AddEntry(signal , "signal", "l");
		tleg->AddEntry(vjets , "v+jets", "l");
		tleg->AddEntry(qcd , "QCD", "l");
	 	tleg->Draw("same");

	 	TText* textChan = doChan(0.12,0.96);
		textChan->Draw();
		TText* textPrelim = doPrelim(0.58,0.96);
		textPrelim->Draw();

	  c1->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_templates.pdf");
	  c1->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_templates.png");
	  delete c1;
	  delete signal;
	  delete vjets;
	  delete qcd;

}

void Fit::normAndColor(TH1D* hist, Sample sample){

	hist->Scale(1./hist->Integral());
	hist->SetMarkerStyle(0);
	hist->SetLineWidth(2);
	hist->SetLineColor(sample.fillColor);
	hist->SetFillColor(kWhite);
}

TH1D* Fit::readHistogram(Sample sample, Variable variable) {

	cout << "plot: " << selection+objName+"/"+variable.name << endl;

	TH1D* plot = (TH1D*) sample.file->Get(selection+objName+"/"+variable.name+"_2btags");
	TH1D* plot2 = (TH1D*) sample.file->Get(selection+objName+"/"+variable.name+"_3btags");
	TH1D* plot3 = (TH1D*) sample.file->Get(selection+objName+"/"+variable.name+"_4orMoreBtags");

	plot->Add(plot2);
	plot->Add(plot3);

	plot->SetFillColor(sample.fillColor);
	plot->SetLineColor(sample.lineColor);

	if(Globals::addOverFlow)
		addOverFlow(plot, variable);

	plot->Rebin(variable.rebinFact);

	return plot;
}

void Fit::addOverFlow(TH1D* overflow, Variable variable){

	if(variable.minX > -0.1){
		int bin = variable.maxX/overflow->GetBinWidth(1);
		double error;

		double overflow_val = overflow->IntegralAndError(bin, overflow->GetNbinsX()+1, error);

		overflow->SetBinContent(bin, overflow_val);
		overflow->SetBinError(bin, error);
	}
}


void Fit::readHistos(AllSamples samples, Variable variable){

	TH1D* data = readHistogram(*samples.single_mu_data, variable);
	TH1D* ttbar = readHistogram(*samples.ttbar, variable);
	TH1D* single_t = readHistogram(*samples.single_t, variable);
	TH1D* signal = (TH1D*) ttbar->Clone("signal");
	signal->Add(single_t);
	TH1D* vjets = readHistogram(*samples.vjets, variable);
	TH1D* qcd = readHistogram(*samples.qcd, variable);


	samples.single_mu_data->SetHisto(data);
	samples.ttbar->SetHisto(ttbar);
	samples.single_t->SetHisto(single_t);
	samples.signal->SetHisto(signal);
	samples.vjets->SetHisto(vjets);
	samples.qcd->SetHisto(qcd);
}

THStack* Fit::buildStack(AllSamples samples, Variable variable){

	THStack *hs = new THStack("hs","test");

	hs->Add(samples.qcd->histo);
	hs->Add(samples.vjets->histo);
	hs->Add(samples.single_t->histo);
	hs->Add(samples.ttbar->histo);

	return hs;
}

TH1D* Fit::allMChisto(AllSamples samples, Variable variable){

	TH1D *allMC = (TH1D*)samples.ttbar->histo->Clone("ratio plot");

	allMC->Add(samples.qcd->histo);
	allMC->Add(samples.vjets->histo);
	allMC->Add(samples.single_t->histo);

	return allMC;
}

void Fit::standardPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable){
	//Style
	TdrStyle style;
	style.setTDRStyle();

	//draw histos to files
	TCanvas *c1 = new TCanvas("Plot","Plot",900, 600);

	data->Draw();
	hs->Draw("hist");

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
		c1->SaveAs("Plots/"+folder+"/"+objName+"/Log/"+variable.name+".png");
		c1->SaveAs("Plots/"+folder+"/"+objName+"/Log/"+variable.name+".pdf");
	}else{
		c1->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+".png");
		c1->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+".pdf");
	}

	delete c1;
	delete leg;
	delete textChan;
	delete textPrelim;
}

void Fit::ratioPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable){
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

//	Will need to see if this works in other situations
	if(variable.minX < 0. && ratio->GetBinWidth(1)){
		ratio->SetAxisRange(variable.minX+ratio->GetBinWidth(1)/2, fabs(variable.minX+ratio->GetBinWidth(1)/2));
	}else{
		ratio->SetAxisRange(0, variable.maxX-ratio->GetBinWidth(1)/2);
	}

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
		c2->SaveAs("Plots/"+folder+"/"+objName+"/Log/"+variable.name+"_ratio.png");
		c2->SaveAs("Plots/"+folder+"/"+objName+"/Log/"+variable.name+"_ratio.pdf");
	}else{
		c2->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_ratio.png");
		c2->SaveAs("Plots/"+folder+"/"+objName+"/"+variable.name+"_ratio.pdf");
	}

	delete c2;
	delete leg;
	delete textChan;
	delete textPrelim;
}

TH1D* Fit::hashErrors(AllSamples samples, Variable variable){
	TH1D * hashErrors = allMChisto(samples, variable);

	hashErrors->SetFillColor(kBlack);
	hashErrors->SetFillStyle(3354);
	hashErrors->SetMarkerSize(0.);
	hashErrors->SetStats(0);

	return hashErrors;
}

TLegend* Fit::legend(AllSamples samples){

		TLegend *tleg;
		tleg = new TLegend(0.75,0.75,0.85,0.9);
		tleg->SetTextSize(0.04);
		tleg->SetBorderSize(0);
		tleg->SetFillColor(10);
		tleg->AddEntry(samples.single_mu_data->histo , "2012 data", "lpe");
		tleg->AddEntry(samples.ttbar->histo , "t#bar{t}", "f");
//		tleg->AddEntry(samples.single_t->histo, "single top"      , "f");
		tleg->AddEntry(samples.vjets->histo , "v+jets", "f");
		tleg->AddEntry(samples.qcd->histo, "QCD"      , "f");

		return tleg;
}

TText* Fit::doChan(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "#mu, #geq 4 jets, #geq 2 btags";

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

TText* Fit::doPrelim(double x_pos,double y_pos){

	  ostringstream stream;
	  stream  << "CMS Preliminary, L = "+Globals::lumi;

	  TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
	  text->SetNDC(true);
	  text->SetTextFont(62);
	  text->SetTextSize(0.045);  // for thesis

	  return text;
}

void Fit::setSelection(TString sel_name){
	selection = sel_name;
}

} /* namespace std */
