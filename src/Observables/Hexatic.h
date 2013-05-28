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

#ifndef HEXATIC_H
#define HEXATIC_H

#include "ObservableBase.h"

complex<double>
Conj21(complex<double> a);

class HexaticClass : public ObservableClass
{
protected:
  int q;
  bool Centroid;
  double DistCutoff;
  Array<complex<double>,1> ParticleOrder;
  void ReadGrid(IOSectionClass &in);
  complex<double> OrderParamater(int slice,int ptcl);
  complex<double> OrderParamater(Array<dVec,1> centroidPos,int ptcl);
  void CalculateCentroid_parallel();
  void Accumulate_old();
  LinearGrid grid;
  Array<complex<double>,1> Histogram;
  Array<double,1> HistDouble;
  Array<double,1> HistSum;
  Array<dVec,1>  CentroidPos;
  ObservableVecDouble1 HexaticRealVar;
  ObservableVecDouble1 HexaticImagVar;
  int NumSamples;
public:
  void Accumulate();
  void WriteInfo();
  void WriteBlock();
  void Read(IOSectionClass &in);
  HexaticClass(PathDataClass &pathData, IOSectionClass &ioSection) :
    ObservableClass (pathData, ioSection),
    HexaticRealVar("HexaticReal",IOSection,pathData.Path.Communicator),
    HexaticImagVar("HexaticImag",IOSection,pathData.Path.Communicator),
    q(6), DistCutoff(2.56)
    {}

};

#endif
