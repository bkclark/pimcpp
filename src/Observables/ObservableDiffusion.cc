/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
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
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "ObservableDiffusion.h"
#include "../Moves/MoveUtils.h"

//dVec v;

////////////////////////////////////////
///      ObsDiffusionClass           ///
////////////////////////////////////////

void ObsDiffusionClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("dumpFrequency",dumpFrequency));
  //in.CloseSection();
}

void ObsDiffusionClass::WriteInfo()
{
  ObservableClass::WriteInfo();
  //IOSection.NewSection("grid");
  //grid.Write(IOSection);
  //IOSection.CloseSection();

  //int numBins = grid.NumPoints-1;
  //Array<double,1> r(numBins);
  //for (int i=0; i<numBins; i++) {
  //  double ra = grid(i);
  //  double rb = grid(i+1);
  //  r(i) = 0.75 * (rb*rb*rb*rb-ra*ra*ra*ra)/(rb*rb*rb-ra*ra*ra);
  //}
  //IOSection.WriteVar("x", r);
  IOSection.WriteVar("xlabel", "t");
  IOSection.WriteVar("ylabel", "MSD");
  IOSection.WriteVar("Type","Autocorrelation");
}

void ObsDiffusionClass::WriteBlock()
{
}

void ObsDiffusionClass::LocalWriteBlock()
{
	if (Path.Communicator.MyProc()==0){
  	if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
	}
	HistVar.Write(TotalMSD);
	TimeStepVar.Write(TotalCounts);
}


void ObsDiffusionClass::Print()
{
/*
  for (int i=0; i<(grid.NumPoints-1); i++)
    {
      double r1 = grid(i);
      double r2 = grid(i+1);
      double r = 0.5*(r1+r2);
      double vol = 4.0*M_PI/3 * (r2*r2*r2-r1*r1*r1);
      double gofr = Histogram(i)/(TotalCounts - WaitToFill);// / (vol*TotalCounts);
      fprintf (stderr, "%1.12e %1.12e\n", r, gofr);
    }
		*/
}

void ObsDiffusionClass::Accumulate()
{
	double norm = 1.0/(totalMol*totalSlices);
  TotalCounts++;
	TotalMSD = 0.0;
  // loop over slices
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    /// loop over molecules 
    for (int mol=0;mol<PathData.Mol.NumMol();mol++){
			dVec rt = PathData.Path(slice,mol) - R0(slice,mol);
			double RSqMag = dot(rt, rt);
			TotalMSD += RSqMag;
    }
  }
	TotalMSD *= norm;
	if(TotalCounts%dumpFrequency == 0){
  	LocalWriteBlock();
	}
}

void ObsDiffusionClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
	totalSlices = PathData.NumTimeSlices()-1;
	totalMol = PathData.Mol.NumMol();
	R0.resize(totalMol, totalSlices);
  for (int slice=0;slice<totalSlices;slice++) {
    /// loop over molecules 
    for (int mol=0;mol<totalMol;mol++){
			dVec r0 = PathData.Path(slice,mol);
			R0(mol,slice) = r0;
		}
	}
}
