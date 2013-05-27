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

#ifndef PIMC_CLASS_H
#define PIMC_CLASS_H


#include "PathDataClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"
#include "LoopClass.h"
#include "RunInfoClass.h"
#include "QMCWrapper.h"

class PIMCClass
{

public:
  std::list<MoveClass*> Moves;
  std::list<ObservableClass*> Observables;
  void ReadMoves(IOSectionClass &in);
  void ReadObservables(IOSectionClass &in);
  void ReadAlgorithm(IOSectionClass &in);
  void WriteSystemInfo();
  void CreateOutFile(IOSectionClass &in);
  string OutFileName;
  IOSectionClass OutFile;
  LoopClass Algorithm;
  RunInfoClass RunInfo;
  // QMC Wrapper
  QMCWrapperClass* QMCWrapper;
public:
  PathDataClass PathData;
  bool Read(IOSectionClass &in);
  void Run();
  void Dummy();
  PIMCClass() : Algorithm(PathData, OutFile, Moves, Observables)
  {
    RunInfo.ProgramName="pimc++";
  }
};


#endif
