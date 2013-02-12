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

#ifndef BISECTION_SPHERE_BLOCK_CLASS_H
#define BISECTION_SPHERE_BLOCK_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "PermuteStage.h"
#include "CoupledPermuteStage.h"
#include "BisectionSphereStage.h"
#include "../Observables/ObservableVar.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionSphereBlockClass : public MultiStageClass
{
private:
  int NumLevels;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  //  void WriteRatio();
  //  ObservableDouble AcceptanceRatioVar;
public:


  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  BisectionSphereBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out)

  { 
    // do nothing for now
  }
};


#endif
