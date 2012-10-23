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


#include "CoupledPermuteStage.h"


double CoupledPermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}


void CoupledPermuteStageClass::InitBlock()
{


}


void CoupledPermuteStageClass::Read (IOSectionClass &in)
{
  double epsilon;
  Array<double,1> gamma;
  assert (in.ReadVar("Gamma", gamma));
  assert (gamma.size()==4);
  assert(in.ReadVar("epsilon",epsilon));
  for (int i=0; i<4; i++) {
    Table1.Gamma[i] = gamma(i);
    Table2.Gamma[i] = gamma(i);
  }
  Table1.epsilon=epsilon;
  Table2.epsilon=epsilon;
}

void CoupledPermuteStageClass::OnlyOdd()
{
  Forw->OnlyOdd=true;
  Rev->OnlyOdd=true;
  Forw->OnlyEven=false;
  Rev->OnlyEven=false;
}

void CoupledPermuteStageClass::OnlyEven()
{
  Forw->OnlyOdd=false;
  Rev->OnlyOdd=false;
  Forw->OnlyEven=true;
  Rev->OnlyEven=true;
}


bool CoupledPermuteStageClass::Attempt (int &slice1, int &slice2, 
				      Array<int,1> &activeParticles, double &prevActionChange)
{

  if (activeParticles(0)==-1){
    Array<int,1> coupledWeight(1);
    Array<int,1> coupledWeightSend(1);
    coupledWeightSend(0)=PathData.Path.SignWeight;
    int myProc=PathData.InterComm.MyProc();
    int numProcs=PathData.InterComm.NumProcs();
    int sendProc=(myProc+1) % numProcs;
    int recvProc=((myProc-1) + numProcs) % numProcs;
    PathData.InterComm.SendReceive (sendProc, coupledWeightSend,
				    recvProc, coupledWeight);
    Array<int,1> doEven(1);
    Array<int,1> dummy(1);
    double myRand=PathData.Path.Random.Local();
    if (myRand>0.5)
      doEven=1;
    else 
      doEven=0;
    if (coupledWeightSend(0)*coupledWeight(0)==1){//Same sign
      if (myProc==0){
	PathData.InterComm.SendReceive(sendProc,doEven,recvProc,dummy);
	if (doEven(0)==1)
	  OnlyOdd();
	else 
	  OnlyEven();
      }
      else if (myProc==1){
	PathData.InterComm.SendReceive(sendProc,dummy,recvProc,doEven);
	if (doEven(0)==1)
	  OnlyEven();
	else 
	  OnlyOdd();
      }
    }
    else {//different sign
      if (myProc==0){
	PathData.InterComm.SendReceive(sendProc,doEven,recvProc,dummy);
	if (doEven(0)==1)
	  OnlyOdd();
	else 
	  OnlyEven();
      }
      if (myProc==1){
	PathData.InterComm.SendReceive(sendProc,dummy,recvProc,doEven);
	if (doEven(0)==1)
	  OnlyOdd();
	else 
	  OnlyEven();
      }
    }
    Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
    int NumPerms = 0;
    // Choose a permutation cycle
    double forwT = Forw->AttemptPermutation();
    activeParticles.resize (Forw->CurrentCycle.Length);
    activeParticles = Forw->CurrentParticles();
    double revT = Rev->CalcReverseProb(*Forw);
    double Tratio = forwT/revT;
    int len=Forw->CurrentCycle.Length;
    if (len % 2 ==0){
      PathData.Path.SignWeight=PathData.Path.SignWeight*-1;
    }
    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    double psi = PathData.Path.Random.Local();

    Array<int,1> currentParticles=Forw->CurrentParticles();  
    double pi_ratio = exp(-actionChange+prevActionChange);
    double acceptProb = min(1.0, pi_ratio/Tratio);
    prevActionChange = actionChange;
    return (acceptProb > psi);
  }
  else
    return true;    
}
