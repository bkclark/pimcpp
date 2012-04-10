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

#ifndef SPECIFIC_HEAT_H
#define SPECIFIC_HEAT_H

#include "ObservableBase.h"

class SpecificHeatClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum;//, TotalActionSum, ExpTotalActionSum, TIP5PSum;
  double WeightA;
  double EA;
  double EB;
  double WeightB;
  ObservableDouble WeightAVar,WeightBVar,EAVar,EBVar;
  //  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
  //    VShortVar, VLongVar;//, TotalActionVar, ExpTotalActionVar, TIP5PVar;
  //ObservableDouble TotalWeightVar;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  vector<ActionBaseClass*> OtherActions;
  vector<ObservableDouble*> OtherVars;
  vector<double> OtherSums;
  int numEnergies;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  SpecificHeatClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      WeightAVar("WeightA",IOSection,myPathData.Path.Communicator),
      WeightBVar("WeightB",IOSection,myPathData.Path.Communicator),
      EAVar("EA",IOSection,myPathData.Path.Communicator),
      EBVar("EB",IOSection,myPathData.Path.Communicator)
    // TotalActionVar ("TotalAction",IOSection,myPathData.Path.Communicator),
    // ExpTotalActionVar ("ExpTotalAction",IOSection,myPathData.Path.Communicator)
    // TIP5PVar  ("TIP5P",IOSection,myPathData.Path.Communicator)
  {
    EA=0.0;
    EB=0.0;
    WeightA=0.0;
    WeightB=0.0;
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
//     TotalActionSum = 0.0;
//     ExpTotalActionSum = 0.0;
    NumSamples = 0;
    TimesCalled=0;
    //    TIP5PSum = 0;
		OtherActions.resize(0);
		OtherVars.resize(0);
		OtherSums.resize(0);
  }
};

#endif 
