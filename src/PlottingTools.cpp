//============================================================================
// Name        : PlottingTools.cpp
// Author      : P Symonds
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "../interface/Objects/MET.h"

using namespace std;

int main() {
	//Have to make a dummy histogram so that everything works. Silly root
	TH1F* dummy = new TH1F("dummy", "dummy", 10, 0, 1);
	delete dummy;

	AllSamples samples;

	MET met;
	met.allPlots(samples);

	return 0;

}
