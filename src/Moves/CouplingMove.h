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

#ifndef COUPLINGMOVE_CLASS_H
#define COUPLINGMOVE_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "MultiStage.h"




/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class CouplingMoveClass : public MultiStageClass
{
 private:
  //  EndType Open;
  int SpeciesNum;
 public:
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass &moveInput);
  void WriteRatio();
  void MakeMove();
  CouplingMoveClass(PathDataClass &myPathData, IOSectionClass iosection) :
    MultiStageClass(myPathData,iosection) 
  { 
    /* Do nothing for now. */ 
  }
};


#endif
