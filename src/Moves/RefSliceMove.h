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

#ifndef REF_SLICE_MOVE_H
#define REF_SLICE_MOVE_H

#include "../PathDataClass.h"
#include "MoveBase.h"
#include "PermuteStage.h"
#include "BisectionStage.h"
#include "DisplaceMove.h"

/// This move, inherited from ParticleMoveClass, performs a set of
/// bisection stages over a set of time slices which contains the
/// reference slice.  The processor which owns the reference slice is
/// the only one that actually moves particles, while the other
/// processors simply compute the change in the nodal action.
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class RefSliceMoveClass : public MultiStageClass
{
private:
  /// Number of bisection stage levels
  int NumLevels;

  int NodeAccept, NodeReject;

  /// Holds the current master processor
  int MasterProc;

  /// The species this particular instance is working on
  int SpeciesNum;
  StageClass *PermuteStage;

  /// This function checks to see if we should accept based on the
  /// change in the node action
  bool NodeCheck();

  /// This function is called if I have the reference slice
  void MakeMoveMaster();

  /// This move is called if I don't have the reference slice
  void MakeMoveSlave();
  bool DoSlaveMoves;

public:
  /// Read in the parameters this class needs from the input file.
  void Read(IOSectionClass &in);
  void WriteRatio();  

  /// Override base class MakeMove to do a block of moves
  void MakeMove();
  inline double AcceptanceRatio() 
  {
    return (double)NodeAccept/(double)(NodeAccept+NodeReject);
  }

  RefSliceMoveClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out)
  { 
    NodeAccept = NodeReject = 0;
    // do nothing for now
  }
};



#endif
