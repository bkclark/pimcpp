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

#ifndef TIMEHEXATIC_H
#define TIMEHEXATIC_H

#include <map>
#include "ObservableBase.h"


class TimeHexaticClass : public ObservableClass
{
protected:
  Array< complex<double> ,2> OPArray;
  Array<dVec , 2> PosArray;
  //  Array<map<int,dVec > >,2>  uj_minus_ujp; // time x ptcl (then list of ptcls you care about w/ the result)
  
  int q;
  bool Centroid;
  double DistCutoff;

  void ReadGrid(IOSectionClass &in);
  void BuildDistanceTable(int time1, int time2);
  ///  complex<double> OrderParamater(int slice,int ptcl);
  complex<double> OrderParamater(Array<dVec,1> &centroidPos,int ptcl);
  ///  void Accumulate_old();
  LinearGrid grid;
  void ProduceTimeMatrix(int slice);
  Array<double,2> DistanceTable;
  
  ObservableVecDouble2 CentroidPosVar;
  ObservableVecDouble2 gofrDispVar;
  ObservableVecDouble2 TimeDispVar;
  ObservableVecDouble1 NumStepVar;
  Array<complex<double> ,2> TimeDisp; // r x time
  Array<double ,2> gofrDisp; // r x time 
  Array<double ,2> TimeDisp_double;
  int TotalCurrentData;
  bool FullyWrapped;
  int NumSamples;
  Array<double,1> NumStepArray;
  int TimeResolution;
  int TimeWindow;
  int NumStoredMCSteps;
  int CurrTime;
  pair<complex<double> ,double>  CalcValue(int time1, int time2,int ptcl1, int ptcl2);
  int numGridPoints;
  double unit2angle(double x,double y);


public:
  void Accumulate();
  Array<dVec,1> CentroidPos;
  Array<double,2> CentroidPos_write;
  Array<complex<double>,1> ParticleOrder;
  void CalculateCentroid();
  void CalculateCentroid_parallel();
  void WriteInfo();
  void WriteBlock();
  void Read(IOSectionClass &in);
  TimeHexaticClass(PathDataClass &pathData, IOSectionClass &ioSection) :
    ObservableClass (pathData, ioSection), 
    TimeDispVar("TimeDisp",IOSection,pathData.Path.Communicator),
    gofrDispVar("gofr",IOSection,pathData.Path.Communicator),
    CentroidPosVar("CentroidPos",IOSection,pathData.Path.Communicator),
    NumStepVar("NumStep",IOSection,pathData.Path.Communicator),
    DistCutoff(2.56), CurrTime(0),TotalCurrentData(-1),FullyWrapped(false),
    q(6)
    {

    }

};

#endif
