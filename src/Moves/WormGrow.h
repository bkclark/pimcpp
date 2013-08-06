/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#ifndef WORM_GROW_MOVE_H
#define WORM_GROW_MOVE_H

#include "../PathDataClass.h"
#include "WormStage.h"
#include "MultiStage.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormGrowMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  WormStageClass WormGrowStage;

	// stuff for debugging
	ofstream out;
	Array<double,2> wormBin;
	int numBins, writeFreq, count;
	double dr;

public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  int MaxGrowth;
  WormGrowMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    WormGrowStage(pathData, outSection)
  {
    DumpFreq = 20;

		// all for debugging
		count = 0;
		writeFreq = 10;
		numBins = 200;
		dr = 0.1;
		wormBin.resize(numBins,2);
		wormBin = 0;
		for(int w=0; w<numBins; w++)
			wormBin(w,0) = dr*w;

	}
};


#endif
