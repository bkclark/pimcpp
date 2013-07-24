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

#ifndef TIMELINDENMAN_H
#define TIMELINDENMAN_H

#include <map>
#include "ObservableBase.h"


class TimeLindenmanClass : public ObservableClass
{
protected:
  Array< map<int,dVec>,2> uj_minus_ujp;
  //  Array<map<int,dVec > >,2>  uj_minus_ujp; // time x ptcl (then list of ptcls you care about w/ the result)
  
  
  bool Centroid;
  double DistCutoff;
  Array<complex<double>,1> ParticleOrder;
  void ReadGrid(IOSectionClass &in);
  ///  complex<double> OrderParamater(int slice,int ptcl);
  ///  complex<double> OrderParamater(Array<dVec,1> centroidPos,int ptcl);
  ///  void Accumulate_old();
  LinearGrid grid;
  Array<complex<double>,1> Histogram;
  Array<double,1> HistDouble;
  Array<double,1> HistSum;
  void ProduceTimeMatrix(int slice);
  ObservableVecDouble1 TimeDispVar;
  ObservableVecDouble1 NumStepVar;
  Array<double,1> TimeDisp;
  int TotalCurrentData;
  bool FullyWrapped;
  int NumSamples;
  Array<double,1> NumStepArray;
  int TimeResolution;
  int TimeWindow;
  int NumStoredMCSteps;
  int CurrTime;
  double CalcValue(int time1, int time2);



public:
  void Accumulate();
  Array<dVec,1> CentroidPos;
  void CalculateCentroid();
  void CalculateCentroid_parallel();
  void WriteInfo();
  void WriteBlock();
  void Read(IOSectionClass &in);
  TimeLindenmanClass(PathDataClass &pathData, IOSectionClass &ioSection) :
    ObservableClass (pathData, ioSection), 
    TimeDispVar("TimeDisp",IOSection,pathData.Path.Communicator),
    NumStepVar("NumStep",IOSection,pathData.Path.Communicator),
    DistCutoff(2.56), CurrTime(0),TotalCurrentData(-1),FullyWrapped(false)
    {

    }

};

#endif
