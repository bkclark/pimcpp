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


#include "PermuteStage.h"


double TablePermuteStageClass::Sample(int &slice1, int &slice2, Array<int,1> &changedParticles)
{
  //do nothing for now
}


void TablePermuteStageClass::Accept()
{
  if (!HaveBeenAcceptedOrRejected){
    int myLen=Forw->CurrentCycle.Length;
    NumAccepted(myLen-1)++;
    NumAttempted(myLen-1)++;
    ///switched for temporary test on cell method ???
    PermuteTableClass* temp;
    if (myLen!=1){
      Forw->UpdateHTable(Forw->CurrentCycle);
      temp=Forw;
      Forw=Rev;
      Rev=temp;
    }

    NeedToRebuildTable = false;
    HaveBeenAcceptedOrRejected = true;
  }
}


void TablePermuteStageClass::Reject()
{
  if (!HaveBeenAcceptedOrRejected){
    int myLen=Forw->CurrentCycle.Length;
    NumAttempted(myLen-1)++;
    if (myLen!=1 && NeedToUpdateHTableOnReject){
      ModeType currMode = GetMode();
      SetMode(OLDMODE);
      CycleClass revCycle;
      revCycle.Length=Forw->CurrentCycle.Length;
      for (int i=0;i<revCycle.Length;i++)
        revCycle.CycleRep(i) = Forw->CurrentCycle.CycleRep(revCycle.Length-(i+1));
      Rev->UpdateHTable(revCycle);
      SetMode(currMode);

    }

    NeedToRebuildTable = false;
    HaveBeenAcceptedOrRejected = true;
  }
}


void TablePermuteStageClass::WriteRatio()
{
   Array<int,1> numAttemptTotal(4);
   Array<int,1> numAcceptTotal(4);
   Array<double,1> ratioTotal(4);
   int totalAttempts=0;
   PathData.Path.Communicator.Sum(NumAttempted,numAttemptTotal);
   PathData.Path.Communicator.Sum(NumAccepted,numAcceptTotal);
   for (int len=0;len<4;len++){
     totalAttempts=totalAttempts+numAcceptTotal(len);
     if (numAttemptTotal(len)!=0)
       ratioTotal(len)=(double)numAcceptTotal(len)/((double)numAttemptTotal(len));
     else
       ratioTotal(len)=0.0;
   }

   if (totalAttempts!=0){
     AcceptanceRatioVar.Write(ratioTotal);
     AcceptanceTotalVar.Write(numAttemptTotal);
     NumAttempted=0;
     NumAccepted=0;
   }

}

void TablePermuteStageClass::InitBlock(int &slice1,int &slice2)
{
  Forw->Slice1 = slice1;
  Forw->Slice2 = slice2;
  Forw->SpeciesNum = SpeciesNum;
  Rev->Slice1 = slice1;
  Rev->Slice2 = slice2;
  Rev->SpeciesNum = SpeciesNum;

  Forw->ConstructHTable();
  Rev->HTable.resize(Forw->HTable.extent(0),Forw->HTable.extent(1));
  for (int i=0; i<Forw->HTable.extent(0); i++) {
    for (int j=0; j<Forw->HTable.extent(1); j++) {
      Rev->HTable(i,j) = Forw->HTable(i,j);
    }
  }
  //  Rev->ConstructHTable();
  NeedToRebuildTable = true;

}


void TablePermuteStageClass::Read (IOSectionClass &in)
{
  Array<double,1> gamma;
  double epsilon;
  assert (in.ReadVar("Gamma", gamma));
  assert (gamma.size()==4);
  assert (in.ReadVar("epsilon",epsilon));
  if (!in.ReadVar("zFocus",zFocus))
    zFocus=false;
  Table1.zfocus=zFocus;
  Table2.zfocus=zFocus;
  double zalpha;
  if (zFocus){
    in.ReadVar("zalpha",zalpha);
    Table1.zalpha=zalpha;
    Table2.zalpha=zalpha;
  }
  for (int i=0; i<4; i++) {
    Table1.Gamma[i] = gamma(i);
    Table2.Gamma[i] = gamma(i);
  }
  Table1.epsilon=epsilon;
  Table2.epsilon=epsilon;

}

//You must call accept or reject after having called attempt once
bool TablePermuteStageClass::Attempt (int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  HaveBeenAcceptedOrRejected = false;
  NeedToUpdateHTableOnReject = false;
  if (activeParticles(0) == -1) { // Initial permutation stage

    // Decide whether or not to rebuild table
    if (PathData.Path.OpenPaths && (PathData.Path.OpenPaths && slice1<=PathData.Path.OpenLink && PathData.Path.OpenLink<=slice2) ||
       (PathData.Path.OpenLink==PathData.Path.NumTimeSlices()-1 && (slice1==0 || slice2==0))) {
      if (NeedToRebuildTable)
        Forw->ConstructCycleTable(SpeciesNum, slice1, slice2, PathData.Path.OpenPtcl);
    } else {
      if (NeedToRebuildTable)
        Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
    }

    // Choose a permutation cycle
    int NumPerms = 0;
    forwT = Forw->AttemptPermutation();
    int len = Forw->CurrentCycle.Length;
    activeParticles.resize(len);
    activeParticles = Forw->CurrentParticles();
    prevActionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    return true;

  } else { // Final permutation stage
    if (Forw->CurrentCycle.Length != 1) {
      NeedToUpdateHTableOnReject = true;
      Rev->UpdateHTable(Forw->CurrentCycle);
    }
    double revT = Rev->CalcReverseProb(*Forw);

    int len = Forw->CurrentCycle.Length;
    if (len % 2 == 0)
      PathData.Path.SignWeight = PathData.Path.SignWeight*-1;

    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    Array<int,1> currentParticles=Forw->CurrentParticles();  
    //    double pi_ratio = exp(-actionChange+prevActionChange);
    double pi_ratio = exp(-actionChange);
    double Tratio = forwT/revT;
    double acceptProb = min(1.0, pi_ratio/Tratio);
    prevActionChange = actionChange;

    double psi = PathData.Path.Random.Local();
    bool toAccept=(acceptProb > psi);
    return toAccept;

  }
}
