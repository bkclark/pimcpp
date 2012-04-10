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

#ifndef CORRELATED_BISECTION_BLOCK_CLASS_H
#define CORRELATED_BISECTION_BLOCK_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "PermuteStage.h"
#include "CoupledPermuteStage.h"
#include "BisectionStage.h"
#include "../Observables/ObservableVar.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class CorrelatedBisectionBlockClass : public MultiStageClass
{

private:
  int StepNum;
  int NumLevels;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  PermuteStageClass* PermuteStage;
  Array<BisectionStageClass*,1> BisectionStages;
  list<ActionBaseClass*> BosonActions;
  list<ActionBaseClass*> NodalActions;
  
  FILE *EAout, *EBout, *SAout, *SBout;
  ObservableDouble wAEAvar, wBEBvar, wAvar, wBvar;
  Array<double,1> AllSumVecIn, AllSumVecOut;
  //  ObservableDouble AcceptanceRatioVar;
public:
  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);

  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  CorrelatedBisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out),StepNum(0), 
    wAEAvar("wAEA", IOSection, pathData.Path.Communicator),
    wBEBvar("wBEB", IOSection, pathData.Path.Communicator),
    wAvar  ("wA",   IOSection, pathData.Path.Communicator),
    wBvar  ("wB",   IOSection, pathData.Path.Communicator),
    AllSumVecIn(2), AllSumVecOut(2)
  { 
    EAout = fopen ("EA.dat", "w");
    SAout = fopen ("SA.dat", "w");
    EBout = fopen ("EB.dat", "w");
    SBout = fopen ("SB.dat", "w");
    // do nothing for now
  }
};


#endif
