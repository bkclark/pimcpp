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

#ifndef BISECTION_BLOCK_CLASS_H
#define BISECTION_BLOCK_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "PermuteStage.h"
#include "CoupledPermuteStage.h"
#include "../Observables/ObservableVar.h"
#include "BisectionStage.h"
#include "BisectionJosephsonStage.h"


// #include "../Actions/FixedPhaseActionBCSClass.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionBlockClass : public MultiStageClass
{
private:
  int StepNum;
  int NumLevels, LowestLevel;
  int StepsPerBlock;
  bool HaveRefslice;
  bool Josephson;
  int SpeciesNum;
  int NumAttempted;
  bool UseApproximateHigherLevelAction;
  void ChooseTimeSlices();
  /// If we do not bisect down to the lowest level, interpolate the
  /// paths in imaginary time.
  void MakeStraightPaths();
  void PrintTimeSpent();
  void WriteRatio();
  StageClass* PermuteStage;
  //  ObservableDouble AcceptanceRatioVar;
public:
  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  void Read_new(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  void MakeMove();
  inline double AcceptanceRatio() 
  {
    return (double)(NumAccepted)/(double)NumAttempted;
  }
  double TimeSpent;
  BisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    //    PermuteRatioVar("PermuteRatioVar",IOSection,pathData.Path.Communicator),
    MultiStageClass(pathData, out),StepNum(0),TimeSpent(0.0)

  { 
    NumAttempted=0;
    UseApproximateHigherLevelAction=false;
    // do nothing for now
  }
  //  FixedPhaseActionBCSClass FP;
};


#endif
