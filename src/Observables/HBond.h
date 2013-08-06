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

#ifndef HBOND_H
#define HBOND_H

#include "ObservableBase.h"

class HbondClass : public ObservableClass
{
	private:

  LinearGrid grid;
	Array<int,3> AgeTable;
	Array<int,2> BondCount;
	Array<int,1> Histogram;
	Array<int,1> LifetimeHist;
	Array<int,1> AgeHist;
	Array<int,2> Protons;

  double TotalSum;

  ObservableVecDouble1 LifetimeVar;
  ObservableVecDouble1 AgeVar;
  ObservableDouble TotalVar;

  int NumSamples;
  int TimesCalled;
	int TotalSamples;
  int Frequency;
  int dumpFrequency;
	double gridStart, gridEnd;
	int numGridPoints;
	int totalSlices;
	int totalMol;
	int overflow;
	int oldest;

	// Here we use HBond criteria r_OO < 3.5 and HOH bond angle > 145 deg, after Artacho (2004)
  double OOlimit;
  double HOHangle;

	public:

  bool CheckPair(int slice, int obond, int ohome, int p);
  bool IsHBond(int slice, int OA, int OB);
  void Accumulate();
  void WriteBlock();
  void Initialize();
  void Read(IOSectionClass& in);
  HbondClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      LifetimeVar  ("Lifetimes",  IOSection,myPathData.Path.Communicator),
      AgeVar  ("Age Distribution",  IOSection,myPathData.Path.Communicator)
  {
    Initialize();
  }

};

#endif 
