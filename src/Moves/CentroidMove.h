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

#ifndef CENTROID_MOVE_H
#define CENTROID_MOVE_H

#include "../PathDataClass.h"
#include "MoveBase.h"
#include "PermuteStage.h"
#include "BisectionStage.h"


/// This stage attempts to displace a list of whole paths.  It should
/// only be used for non-permuting particles.
class CentroidStageClass : public CommonStageClass
{
public:
  void Accept();
  void Reject();

  /// This does the actual displacement of the path.  All processors
  /// within a single close must displace by the same amount.
  double Sample (int &slice1, int &slice2, Array <int,1> &activeParticles);
  CentroidStageClass (PathDataClass &pathData, IOSectionClass &out) :
    CommonStageClass (pathData, out)
  {
    // Do nothing for now.
  }
};


/// This move, inherited from ParticleMoveClass, performs a
/// bisective move, followed by a move of the particles
/// centroid back to the original position
class CentroidMoveClass : public MultiStageClass
{
private:
  /// Number of bisection stage levels
  int NumLevels;

  /// Holds the current master processor
  int MasterProc;

  /// The species this particular instance is working on
  int SpeciesNum;
  StageClass *PermuteStage;

  /// Number of attempted moves
  int NumAttempted;

  /// Choose which slices are bisected
  void ChooseTimeSlices();

  /// The displace stage for the centroid to return
  CentroidStageClass DisplaceStage;
public:
  /// Read in the parameters this class needs from the input file.
  void Read(IOSectionClass &in);
  void WriteRatio();

  /// Override base class MakeMove to do a block of moves
  void MakeMove();
  inline double AcceptanceRatio()
  {
    double accept_ratio=(double)(NumAccepted)/(double)NumAttempted;
    return accept_ratio;
  }

  CentroidMoveClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out),
    DisplaceStage(pathData, out)
  {
    NumAccepted = NumAttempted = 0;
  }
};



#endif
