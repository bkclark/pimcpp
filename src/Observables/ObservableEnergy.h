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

#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H

#include "ObservableBase.h"
#include "HistogramClass.h"
class EnergyClass : public ObservableClass
{

private:
  int GetPermNumber();
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum, dUNonlocalSum,Residual;
  //, TotalActionSum, ExpTotalActionSum, TIP5PSum;
  Array<double,1> EnergyVals;
  ObservableVecDouble1 EnergyValsVar;
  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar, dUNonlocalVar,ResidualVar;
  ObservableVecDouble1 VTailSRVar;
  ObservableVecDouble1 EnergyHistogramVar;
  ObservableDouble VTailLRVar;
  ObservableDouble HistStart;
  ObservableDouble HistEnd;
  ObservableDouble NumPoints;
  //, TotalActionVar, ExpTotalActionVar, TIP5PVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  vector<ActionBaseClass*> OtherActions;
  vector<ObservableDouble*> OtherVars;
  vector<double> OtherSums;
  int numEnergies;
  Array<bool,1> CountedAlready;
  Array<int,1> TotalPerm;
  HistogramClass EnergyHistogram;
  Array<double,1> EnergyHistogramSum;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
    TotalVar     ("Total",  IOSection,myPathData.Path.Communicator),
    KineticVar   ("Kinetic",IOSection,myPathData.Path.Communicator), 
    dUShortVar   ("dUShort",IOSection,myPathData.Path.Communicator), 
    dULongVar    ("dULong", IOSection,myPathData.Path.Communicator), 
    NodeVar      ("Node",   IOSection,myPathData.Path.Communicator), 
    VShortVar    ("VShort",IOSection,myPathData.Path.Communicator), 
    VLongVar     ("VLong",IOSection,myPathData.Path.Communicator),
    dUNonlocalVar("dUNonlocal", IOSection,myPathData.Path.Communicator),
    EnergyValsVar("Energy Vals",IOSection,myPathData.Path.Communicator),
    ResidualVar("Residual Energy",IOSection,myPathData.Path.Communicator),
      VTailSRVar("VTail Short Range",IOSection,myPathData.Path.Communicator),
      VTailLRVar("VTail Long Range",IOSection,myPathData.Path.Communicator),
      EnergyHistogramVar("Energy Histogram",IOSection,myPathData.Path.Communicator),
      HistStart("HistStart",IOSection,myPathData.Path.Communicator),
      HistEnd("HistEnd",IOSection,myPathData.Path.Communicator),
      NumPoints("NumPoints",IOSection,myPathData.Path.Communicator)
      
// TotalActionVar ("TotalAction",IOSection,myPathData.Path.Communicator),
// ExpTotalActionVar ("ExpTotalAction",IOSection,myPathData.Path.Communicator)
// TIP5PVar  ("TIP5P",IOSection,myPathData.Path.Communicator)
  {
    TotalSum      = 0.0;
    KineticSum    = 0.0;
    dUShortSum    = 0.0;
    dULongSum     = 0.0;
    NodeSum       = 0.0;
    VShortSum     = 0.0;
    VLongSum      = 0.0;
    dUNonlocalSum = 0.0;
    Residual=0.0;
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

class EnergySignClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum, dUNonlocalSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar, dUNonlocalVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);

  EnergySignClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar     ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar   ("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar   ("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar    ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar      ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar    ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar     ("VLong",IOSection,myPathData.Path.Communicator),
      dUNonlocalVar("dUNonlocal", IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
