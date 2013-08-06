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

#ifndef PAIR_CORRELATION_REWEIGHTING_H
#define PAIR_CORRELATION_REWEIGHTING_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class PairCorrelationReweightingClass : public ObservableClass
{
private:
  ObservableVecDouble1 gofrVar;
  /// Stores number of counts in each bin
  //Array<int,1> Histogram;
  // the histogram stores doubles to accommodate reweighting factor
  Array<double,1> Histogram;
  Array<double,1> Correction;
  double delta_S, NormWeight, CorrectWeight;
  // actions for reweighting
  list<ActionBaseClass*> SampledActions;
  list<ActionBaseClass*> ReweightActions;
  Array<int,1> activeParticles;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  /// The species between which I am calculating the pair correlation
  /// function.
  int Species1, Species2;
  /// This grid defines the bins.  Bin 0 is bounded by 0 on the 
  /// bottom and grid(0) on the top.
  LinearGrid grid;
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  PairCorrelationReweightingClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    gofrVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
  PairCorrelationReweightingClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection), 
    Species1(species1), Species2(species2),
    gofrVar("y", IOSection, myPathData.Path.Communicator)
  { Initialize(); }


};

#endif
