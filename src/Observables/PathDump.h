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

#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"


class PathDumpClass : public ObservableClass
{
private:
  ObservableVecDouble3 PathVar;
  ObservableVecInt1 PermVar;
  ObservableInt OpenLinkVar;
  ObservableInt OpenLinkPtclVar;
  ObservableInt RefLinkVar;
  ObservableVecDouble1 TailLocVar;
  // Variables for node dumps.
  ObservableVecDouble3 NodeVarA, NodeVarB;
  ObservableVecDouble3 RhoVar;
  ObservableVecDouble1 WarpPosVar;
  ObservableInt NodePtclVar, NodeSliceVar;
  bool DumpNodes, DumpRho;
  /// If AllClones is true, then all clones will generate path dumps.
  /// Otherwise only clone 0 will.
  bool AllClones;
  int  NodePtcl;
  int  NodeSlice;
  LinearGrid Xgrid, Ygrid, Zgrid;
  void NodeDump();
  void FreeParticleNodeDump();
  void GroundStateNodeDump();
  void FixedPhaseNodeDump();
  void FindWorstBead(int &slice, int &ptcl);
  Array<double,3> Rho;
public:
  int TimesCalled;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  int DumpFreq;
  PathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection) :
    ObservableClass(myPathData, ioSection),
    PathVar ("Path", IOSection, myPathData.Path.Communicator),
    PermVar ("Permutation", IOSection, myPathData.Path.Communicator),
    OpenLinkVar("OpenLinkSlice",IOSection,myPathData.Path.Communicator),
    TailLocVar("TailLocation",IOSection,myPathData.Path.Communicator),
    OpenLinkPtclVar("OpenPtcl",IOSection,myPathData.Path.Communicator),
    RefLinkVar("RefLink",IOSection,myPathData.Path.Communicator),
    NodeVarA("ANodes",IOSection,myPathData.Path.Communicator),
    NodeVarB("BNodes",IOSection,myPathData.Path.Communicator),
    NodePtclVar ("NodePtcl",  IOSection, myPathData.Path.Communicator),
    NodeSliceVar("NodeSlice", IOSection, myPathData.Path.Communicator),
    WarpPosVar("WarpPos", IOSection, myPathData.Path.Communicator),
    RhoVar("Rho", IOSection, myPathData.Path.Communicator),
    AllClones(false)
  { 
    Name="PathDump";
    TimesCalled=0;
  }
};



#endif
