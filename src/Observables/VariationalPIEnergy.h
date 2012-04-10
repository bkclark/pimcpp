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

#ifndef VARIATIONALPI_ENERGY_H
#define VARIATIONALPI_ENERGY_H



#include "ObservableBase.h"




class VariationalPIEnergyClass : public ObservableClass
{

private:
  double Energy;
  Array<double,2> DetMatrix;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  ObservableDouble EnergyVar;
public:
  double rho(int i,int j);
  double DRho(int i, int j);
  int NumImages;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  VariationalPIEnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      EnergyVar("Energy",IOSection,PathData.Path.Communicator)
  {
    Energy=0.0;
    NumSamples = 0;
    TimesCalled=0;
    NumImages=0;
  }
};


#endif 
