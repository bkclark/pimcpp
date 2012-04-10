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

#ifndef VARIATIONAL_DISPLACE_MOVE_H
#define VARIATIONAL_DISPLACE_MOVE_H

#include "../PathDataClass.h"
#include "MultiStage.h"


/// This stage attempts to displace a list of whole paths.  It should
/// only be used for non-permuting particles.
class VariationalDisplaceStageClass : public CommonStageClass
{
public:
  /// This is the width of the gaussian distribution for the
  /// displacement vector.
  bool Attempt(int &slice1, int &slice2, 
	       Array<int,1> &activeParticles,
	       double &prevActionChange);
  double Sigma;
  /// This does the actual displacement of the path.  All processors
  /// within a single close must displace by the same amount.
  double Sample (int &slice1, int &slice2,
		 Array <int,1> &activeParticles);
  VariationalDisplaceStageClass (PathDataClass &pathData,IOSectionClass &outSection) :
    CommonStageClass (pathData,outSection) 
  {
    // Do nothing for now.
  }
};



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class VariationalDisplaceMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  VariationalDisplaceStageClass VariationalDisplaceStage;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  VariationalDisplaceMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    VariationalDisplaceStage(pathData, outSection)
  {
    DumpFreq = 20;
  }
};


#endif
