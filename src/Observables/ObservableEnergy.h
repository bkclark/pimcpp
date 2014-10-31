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

#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H

#include "ObservableBase.h"
#include "HistogramClass.h"


class EnergyClass : public ObservableClass
{

private:
  ObservableDouble TotalVar, VShortVar, VLongVar, VExtVar, HistStart, HistEnd, NumPoints;
  ObservableVecDouble1 vLong_r0_var, vLong_k0_var, duLong_r0_var, duLong_k0_var, EnergyHistogramVar, PermEnergyVar;

  double TotalSum, VShortSum, VLongSum, VExtSum;
  map<string,double> energies;
  map<string,double> actions;
  map<string,double> ESum;
  map<string,ObservableDouble*> EVar;

  vector<double> PermEnergy, SectorCount;
  Array<double,1> EnergyHistogramSum;
  Array<double,1> vLong_k0, vLong_r0, duLong_k0, duLong_r0;

  int numEnergies;

  HistogramClass EnergyHistogram;

  int NumSamples, TimesCalled, Freq, DumpFreq;
  bool CountPerms;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);

  EnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      TotalVar("Total",IOSection,myPathData.Path.Communicator),
      VShortVar("VShort",IOSection,myPathData.Path.Communicator),
      VLongVar("VLong",IOSection,myPathData.Path.Communicator),
      VExtVar("VExt",IOSection,myPathData.Path.Communicator),
      vLong_k0_var("vLong_k0",IOSection,myPathData.Path.Communicator),
      vLong_r0_var("vLong_r0",IOSection,myPathData.Path.Communicator),
      duLong_k0_var("duLong_k0",IOSection,myPathData.Path.Communicator),
      duLong_r0_var("duLong_r0",IOSection,myPathData.Path.Communicator),
      PermEnergyVar("PermEnergy",IOSection,myPathData.Path.Communicator),
      EnergyHistogramVar("Energy Histogram",IOSection,myPathData.Path.Communicator),
      HistStart("HistStart",IOSection,myPathData.Path.Communicator),
      HistEnd("HistEnd",IOSection,myPathData.Path.Communicator),
      NumPoints("NumPoints",IOSection,myPathData.Path.Communicator)
  {
    TotalSum = 0.0;
    VShortSum = 0.0;
    VLongSum = 0.0;
    VExtSum = 0.0;
    NumSamples = 0.0;
    TimesCalled = 0.0;
  }
};


#endif
