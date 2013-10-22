/*
 * MET.h
 *
 *  Created on: Oct 11, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Object.h"

#ifndef MET_H_
#define MET_H_

namespace std {

class MET: public Object {
public:
	MET();
	virtual ~MET();
	void allPlots(AllSamples samples);
//private:
//	TLegend* legend(AllSamples samples);
//	TText* doPrelim(double x_pos,double y_pos);
//	TText* doChan(double x_pos,double y_pos);
//	TH1D* readGe2bHistogram(Sample sample, Variable variable);
//	THStack* buildStack(AllSamples samples, Variable variable);
//	void savePlot(AllSamples samples, Variable variable);
};

} /* namespace std */
#endif /* MET_H_ */
