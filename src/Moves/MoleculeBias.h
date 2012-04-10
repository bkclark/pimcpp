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

#ifndef MOLECULE_BIAS_MOVES_H
#define MOLECULE_BIAS_MOVES_H

#include "MoleculeMoveBase.h"
#include "MoleculeMove.h"
#include "../Actions/MoleculeInteractionsClass.h"

class MoleculeForceBiasMove : public MolMoveClass
{
  int numGen,numProp;
  bool doTrans, doRot;
  double numInProp;
  MoleculeInteractionsClass* MolAction;

 public:
  double Theta;
  double Step;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  MoleculeForceBiasMove(PathDataClass &myPathData,IOSectionClass outSection);
};

dVec CalcLever(dVec axis, dVec coord);

//class ForceBias : public ParticleMoveClass
//{
// public:
//  int numAccepted,numMoves,numGen,numProp;
//  void MakeMove();
//  inline void Read(IOSectionClass &moveInput)
//    {
//      string typeCheck;
//      assert(moveInput.ReadVar("type",typeCheck));
//      assert(typeCheck=="ForceBias");
//      assert(moveInput.ReadVar("name",Name));
//
//    }
//  void AssignPtcl(int mol,Array<int,1>& activeParticles);
//  void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);
//  dVec Normalize(dVec v);
//  double Dot(dVec vec1, dVec vec2);
//  dVec Cross(dVec v1, dVec v2);
//  double Mag(dVec v);
//  dVec Scale(dVec v, double scale);
//  dVec CalcLever(dVec axis, dVec coord);
//  dVec ArbitraryRotate(dVec axis,dVec coord, double phi);
//  ForceBias(PathDataClass &myPathData,IOSectionClass outSection) : 
//    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0),numGen(0),numProp(0)
//    {
//      ActiveParticles.resize(5);
//    }
//  //  double AcceptanceRatio(int numAccepted,int numMoves);
//  void WriteRatio()
//    {
//      //do nothing for now
//    };
//};

#endif
