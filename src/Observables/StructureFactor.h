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

///These work in parallel except for skmax
#ifndef STRUCTURE_FACTOR_H
#define STRUCTURE_FACTOR_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class StructureFactorClass : public ObservableClass
{
  ObservableVecDouble1 SofkVar;
  /// Stores the total number of counts
  Array<dVec,1> Additionalkvecs;
  Array<complex<double>,3> AdditionalRho_k;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  int TotalCounts;
  ObservableDouble SkMaxVar;
  double SkMax;
  dVec MaxkVec;
  /// Stores number of counts in each bin
  Array<double,1> Sk;
  Array<double,1> rho_k_real;
  Array<double,1> rho_k_imag;
  ObservableVecDouble1 rho_k_realVar;
  ObservableVecDouble1 rho_k_imagVar;
  /// The species between which I am calculating the pair correlation
  /// function.
  int Species1, Species2;
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  void Calculate();
  void Clear();
  StructureFactorClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection),
    SkMaxVar("SkMax",IOSection,myPathData.Path.Communicator),
    rho_k_realVar("RhoK_real",IOSection,myPathData.Path.Communicator),
    rho_k_imagVar("RhoK_imag",IOSection,myPathData.Path.Communicator),
    SofkVar("y", IOSection, myPathData.Path.Communicator)
    {
      TimesCalled=0;
    }
  StructureFactorClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection), 
    Species1(species1), Species2(species2),
    SkMaxVar("SkMax",IOSection,myPathData.Path.Communicator),
    SofkVar("y", IOSection, myPathData.Path.Communicator),
    rho_k_realVar("RhoK_real",IOSection,myPathData.Path.Communicator),
    rho_k_imagVar("RhoK_imag",IOSection,myPathData.Path.Communicator)
    { Initialize(); }
};



#endif
