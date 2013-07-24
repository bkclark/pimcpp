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

#ifndef CLUSTER_MOVE_H
#define CLUSTER_MOVE_H

#include "MoveBase.h"
#include "../Actions/MoleculeInteractionsClass.h"

class LocalFlip : public ParticleMoveClass
{
  MoleculeInteractionsClass* MolAction;
 public:
  int numAccepted,numMoves;
  void MakeMove();
  void Read(IOSectionClass &moveInput);
  void AssignPtcl(int mol,Array<int,1>& activeParticles);
  double MolPairAction(int slice,int m,int n);
 // void RotateMol(int slice,int mol,dVec axis,double phi);
  void RotateMol(int slice,int mol,dVec Q);
  dVec ArbitraryRotate(dVec axis,dVec coord, double phi);
  void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);
  void IntegrityCheck(int slice, Array<int,1> activeParticles);
  LocalFlip(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

class GlobalFlip : public ParticleMoveClass
{
  MoleculeInteractionsClass* MolAction;
 public:
  int numAccepted,numMoves;
  void MakeMove();
  void Read(IOSectionClass &moveInput);
  void AssignPtcl(int mol,Array<int,1>& activeParticles);
  double MolPairAction(int slice,int m,int n);
  void RotateMol(int slice,int mol,dVec Q);
  GlobalFlip(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

#endif
