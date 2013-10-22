/*
 * Sample.h
 *
 *  Created on: Oct 14, 2013
 *      Author: philip
 */
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include <iostream>

#ifndef SAMPLE_H_
#define SAMPLE_H_

namespace std {

class Sample {
public:
	Sample();
	Sample(TString sample, Color_t fColor, Color_t lColor);
	virtual ~Sample();
	TFile* file;
	TH1D* histo;
	void SetHisto(TH1D* plot);
	void SetFillColor(Color_t fcolor); //pub or priv?
	Color_t GetFillColor();
	void SetLineColor(Color_t fcolor); //pub or priv?
	Color_t GetLineColor();
	Color_t fillColor;
	Color_t lineColor;
};

} /* namespace std */
#endif /* SAMPLE_H_ */
