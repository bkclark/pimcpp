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

#include "PermuteTableClass.h"


// double PermuteTableClass::Sample (int &slice1, int &slice2,
// 				  Array<double,1> &activeParticles)
// {
//   if (activeParticles(0) == -1) {
//     // We need to choose the permutation now


//     return 1.0;
//   }
//   else {
//     // We need to compute the transition probability now
    
//     // Construct the reverse table;
    
//     // Calculate the reverse probability
    
//     // Return ratio
//   }
  
// }

    
Array<int,1> PermuteTableClass::CurrentParticles()
{
  Array<int,1> tempPtcl(CurrentCycle.Length);
  int firstPtcl=PathData.Species(SpeciesNum).FirstPtcl;
  for (int i=0;i<CurrentCycle.Length;i++){
    tempPtcl(i)=CurrentCycle.CycleRep(i)+firstPtcl;
  }
  return tempPtcl;
}

void PermuteTableClass::UpdateHTable(const CycleClass &perm)
{
 //  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
//   int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
//   int N = lastPtcl-firstPtcl+1;
//   double tempDistance;
//   for (int i=0; i<N; i++) {
//     tempDistance=HTable(i,perm.CycleRep(0));
//     for (int counter=1;counter<perm.Length;counter++){
//       HTable(i,perm.CycleRep(counter-1))=HTable(i,perm.CycleRep(counter));
//     }
//     HTable(i,perm.CycleRep(perm.Length-1))=tempDistance;
//   }    


//   for (int i=0; i<N; i++) {
//     tempDistance=HTable(perm.CycleRep(0),i);
//     for (int counter=1;counter<perm.Length;counter++){
//       HTable(perm.CycleRep(counter-1),i)=HTable(perm.CycleRep(counter),i);
//     }
//     HTable(perm.CycleRep(perm.Length-1),i)=tempDistance;
//   }    


  double logEps=log(epsilon);
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Path.tau * (double) (Slice2-Slice1);
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  int N = lastPtcl-firstPtcl+1;
  for (int i=0; i<N; i++) {
    dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
    PathData.Path.PutInBox(disp_ii);
    double dist_ii = dot (disp_ii,disp_ii);
    for (int jC=0; jC<perm.Length; jC++) {
      int j=perm.CycleRep(jC);
      dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
      PathData.Path.PutInBox(disp_ij);
      double dist_ij = dot (disp_ij,disp_ij);
      double toExp=(-dist_ij + dist_ii)*fourLambdaBetaInv;
      if (toExp<logEps)
        HTable(i,j)=0.0;
      else
        HTable(i,j)=exp(toExp);
    }
  }

  for (int iC=0; iC<perm.Length; iC++) {
    int i=perm.CycleRep(iC);
    dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
    PathData.Path.PutInBox(disp_ii);
    double dist_ii = dot (disp_ii,disp_ii);
    for (int j=0; j<N; j++) {
      dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
      PathData.Path.PutInBox(disp_ij);
      double dist_ij = dot (disp_ij,disp_ij);
      double toExp=(-dist_ij + dist_ii)*fourLambdaBetaInv;
      if (toExp<logEps)
        HTable(i,j)=0.0;
      else
        HTable(i,j)=exp(toExp);
    }
  }

}

void PermuteTableClass::ConstructHTable()
{
  double logEps=log(epsilon);
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Path.tau * (double) (Slice2-Slice1);
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  int N = lastPtcl-firstPtcl+1;
  HTable.resize(N,N);
  //  cerr<<"N: "<<N<<" "<<fourLambdaBetaInv<<endl;
  for (int i=0; i<N; i++) {
    dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
    PathData.Path.PutInBox(disp_ii);
    double dist_ii = dot (disp_ii,disp_ii);
    for (int j=0; j<N; j++) {
      dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
      PathData.Path.PutInBox(disp_ij);
      double dist_ij = dot (disp_ij,disp_ij);
      double toExp=(-dist_ij + dist_ii)*fourLambdaBetaInv;
      if (toExp<logEps)
        HTable(i,j)=0.0;
      else
        HTable(i,j)=exp(toExp);
    // Now we should really sort with respect to j, but we won't for
    // now because we are far too lazy and it's a Friday.
    }
  }
}
///// This constructs the Htable for the reverse move from the Htable
///// for the forward move.  i.e. Constructs Hrev(i,j) = H(P(i),P(j)).
// void PermuteTableClass::PermuteHTable(const CycleClass &myPerm,
// 				      const PermuteTableClass &forwardTable)
// {
//   int N=HTable.extent(0);
//   Array<int,1> P(PathData.NumParticles());
//   myPerm.CanonicalPermRep(P);
//   for (int i=0; i<N; i++) {
//     int P_i = P(i);
//     for (int j=0; j<N; j++) {
//       int P_j = P(j);
//       HTable(P_i,P_j) = forwardTable.HTable(i,j);
//     }
//   }
 
// } 

 
// Actually applies my cyclic permutation to a path at a given time
// slice and the permutation vector in the path
void CycleClass::Apply(PathClass &path, int firstPtcl, int slice)
{
  SetMode(NEWMODE);
  dVec tempPos = path(slice, CycleRep(0)+firstPtcl);
  int tempPtcl = path.Permutation(CycleRep(0)+firstPtcl);
  for(int i=0;i<Length-1;i++) {
    path.SetPos(slice, CycleRep(i)+firstPtcl, path(slice,CycleRep(i+1)+firstPtcl));
    path.Permutation(CycleRep(i)+firstPtcl) = path.Permutation(CycleRep(i+1)+firstPtcl);
  }
  path.SetPos(slice,CycleRep(Length-1)+firstPtcl,tempPos);
  path.Permutation(CycleRep(Length-1)+firstPtcl) = tempPtcl;
}

double PermuteTableClass::AttemptPermutation()
{
  double xi = PathData.Path.Random.Local();
  int index = FindEntry(xi);
  CurrentCycle = CycleTable(index);
  if (CurrentCycle.CycleRep[0]<0 || CurrentCycle.CycleRep[0]>10000){
    PrintTable();
  }
  // Now, apply the permutation to the Path
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  CurrentCycle.Apply(PathData.Path,firstPtcl,Slice2);
  return (CurrentCycle.P * NormInv);
}

double PermuteTableClass::CalcReverseProb(const PermuteTableClass &forwardTable)
{
  // Reverse probabily for a single ptcl move is the same as that for the
  // forward move.
  if (forwardTable.CurrentCycle.Length == 1)
    return (forwardTable.CurrentCycle.P*forwardTable.NormInv);
  //We reconstruct things from scratch to make sure we do the right
  //thing. We can try to incrementally make this faster.
  ConstructCycleTable(forwardTable.SpeciesNum, forwardTable.Slice1, forwardTable.Slice2, forwardTable.ExcludeParticle);
  int len=forwardTable.CurrentCycle.Length;
  return ((Gamma[len-1]*Gamma[len-1])/(forwardTable.CurrentCycle.P)*NormInv);
}


void PermuteTableClass::Read(IOSectionClass &inSection)
{
  Array<double,1> tempGamma;
  assert(inSection.ReadVar("Gamma",tempGamma));
  assert(tempGamma.size() == 4);
  for (int counter=0;counter<tempGamma.size();counter++){
    Gamma(counter)=tempGamma(counter);
  }
  assert(inSection.ReadVar("epsilon",epsilon));
  inSection.ReadVar("zfocus",zfocus);
  if (zfocus)
    assert(inSection.ReadVar("zalpha",zalpha));
  else
    zalpha=0.0;
  //  assert(inSection.ReadVar("SpeciesNum",SpeciesNum));
}


void PermuteTableClass::PrintTable() const
{
  for (int i=0; i<NumEntries; i++) {
    const CycleClass &cycle = CycleTable(i);
    cerr << "Length = " << cycle.Length << endl;
    cerr << "Cycle = [";
    for (int j=0; j<cycle.Length; j++)
      cerr << cycle.CycleRep[j] << " ";
    cerr << "]\n";
    cerr << "P = " << cycle.P << " C = " << cycle.C << endl;
    cerr << endl;
  }
}

void PermuteTableClass::ConstructCycleTable(int speciesNum, int slice1, int slice2)
{
  ConstructCycleTable(speciesNum,slice1,slice2,-1);
}

void PermuteTableClass::ConstructCycleTable(int speciesNum, int slice1, int slice2, int excludeParticle)
{
  /// HACK HACK HACK HACK
   if (PathData.Path.Species(SpeciesNum).GetParticleType()==FERMION)
     if (PathData.Path.UseNodeImportance)
       ConstructBosonCycleTable(speciesNum,slice1,slice2,excludeParticle);
     else
       ConstructFermionCycleTable(speciesNum,slice1,slice2);
   else
     ConstructBosonCycleTable(speciesNum,slice1,slice2,excludeParticle);
}

///the Gamma's must all be greater then 1 or this won't do the correct
///thing.  Allows you to exclude a particle when producing the table
///(say for example the open particle).
void PermuteTableClass::ConstructFermionCycleTable(int speciesNum, int slice1, int slice2)
{
  Slice1=slice1;
  Slice2=slice2;
  SpeciesNum=speciesNum;
  //  ConstructHTable();
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  CycleTable.resize(TableSize);
  NumEntries = 0;
  int N=HTable.extent(0);
  for (int i=0; i<N; i++) {
    ///Single particle move
    double hprod=1.0;
    double hprod2=1.0;
    double hprod3=1.0;
    CycleClass tempPerm;
    tempPerm.Length=1;
    tempPerm.CycleRep[0]=i;
    tempPerm.P=Gamma(0);
    tempPerm.C=tempPerm.P;
    if (NumEntries!=0){
      tempPerm.C+=CycleTable(NumEntries-1).C;
    }
    AddEntry(tempPerm);
    for (int j=i+1; j<N; j++) {// 2 and higher cycles
      //2 cycle permutations
      tempPerm.CycleRep[1]=j;
      hprod=HTable(i,j);
      tempPerm.Length=3;
      for (int k=i+1;k<N;k++){//3 and higher cycles
        //3 cycle permutations
        if (k!=j){
          tempPerm.CycleRep[2]=k;
          hprod2=hprod*HTable(j,k);
          tempPerm.P=Gamma(2)*hprod2*HTable(k,i);
          tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;
          AddEntry(tempPerm);
        }
      }
    }
  }


  Norm = CycleTable(NumEntries-1).C;
  NormInv = 1.0/Norm;

}


///the Gamma's must all be greater then 1 or this won't do the correct
///thing.  Allows you to exclude a particle when producing the table
///(say for example the open particle).
void PermuteTableClass::ConstructBosonCycleTable(int speciesNum, int slice1, int slice2, int particleExclude)
{
  Slice1=slice1;
  Slice2=slice2;
  SpeciesNum=speciesNum;
  ExcludeParticle=particleExclude;
  //  ConstructHTable();
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  //  CycleTable.resize(TableSize);
  NumEntries = 0;
  int N=HTable.extent(0);
  for (int i=0; i<N; i++) {
    ///Single particle move
    if (i!=particleExclude){
      double hprod=1.0;
      double hprod2=1.0;
      double hprod3=1.0;
      CycleClass tempPerm;
      tempPerm.Length=1;
      tempPerm.CycleRep[0]=i;
      tempPerm.P=Gamma(0);
      if (zfocus)
        tempPerm.P=tempPerm.P*exp(-zalpha*PathData.Path(slice1,i+firstPtcl)[2]*PathData.Path(slice1,i+firstPtcl)[2]+-zalpha*PathData.Path(slice2,i+firstPtcl)[2]*PathData.Path(slice2,i+firstPtcl)[2]);
      tempPerm.C=tempPerm.P;
      if (NumEntries!=0){
        tempPerm.C+=CycleTable(NumEntries-1).C;
      }
      AddEntry(tempPerm);
      for (int j=i+1; j<N; j++) {// 2 and higher cycles
      if (j!=particleExclude && i!=particleExclude){
          //2 cycle permutations
          tempPerm.CycleRep[1]=j;
          hprod=HTable(i,j);
          tempPerm.Length=2;
          tempPerm.P=Gamma(1)*hprod*HTable(j,i);
          if (zfocus)
            tempPerm.P=tempPerm.P*exp(-zalpha*PathData.Path(slice1,j+firstPtcl)[2]*PathData.Path(slice1,j+firstPtcl)[2]+-zalpha*PathData.Path(slice2,j+firstPtcl)[2]*PathData.Path(slice2,j+firstPtcl)[2]);
          tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;
          if (tempPerm.P > epsilon)
            AddEntry(tempPerm);
          if (hprod * Gamma(2)*Gamma(3)>epsilon) {
            for (int k=i+1;k<N;k++){//3 and higher cycles
              if (k!=particleExclude){
                //3 cycle permutations
                if (k!=j){
                  tempPerm.CycleRep[2]=k;
                  hprod2=hprod*HTable(j,k);
                  tempPerm.Length=3;
                  tempPerm.P=Gamma(2)*hprod2*HTable(k,i);
                  if (zfocus)
                    tempPerm.P=tempPerm.P*exp(-zalpha*PathData.Path(slice1,k+firstPtcl)[2]*PathData.Path(slice1,k+firstPtcl)[2]+-zalpha*PathData.Path(slice2,k+firstPtcl)[2]*PathData.Path(slice2,k+firstPtcl)[2]);
                  tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;
                  if (tempPerm.P > epsilon)
                    AddEntry(tempPerm);
                  if (hprod2 *Gamma(3) > epsilon) {
                    for (int l=i+1;l<N;l++)
                      if (l!=particleExclude && (l!=j) && (l!=k)){
                        hprod3=hprod2*HTable(k,l);
                        tempPerm.CycleRep[3]=l;
                        tempPerm.Length=4;
                        tempPerm.P=Gamma(3)*hprod3*HTable(l,i);
                          if (zfocus)
                            tempPerm.P=tempPerm.P*exp(-zalpha*PathData.Path(slice1,l+firstPtcl)[2]*PathData.Path(slice1,l+firstPtcl)[2]+-zalpha*PathData.Path(slice2,l+firstPtcl)[2]*PathData.Path(slice2,l+firstPtcl)[2]);
                        tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;
                        if (tempPerm.P > epsilon)
                          AddEntry(tempPerm);
                      }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
//   if (OnlyOdd || PathData.Path.Species(SpeciesNum).GetParticleType()==FERMION){
//     int currSpot=0;
//     for (int counter=0;counter<NumEntries;counter++){
//       if (CycleTable(counter).Length % 2 ==1){
// 	CycleTable(currSpot)=CycleTable(counter);
// 	currSpot++;
//       }
//     }
//     NumEntries=currSpot;
//     CycleTable(0).C=CycleTable(0).P;
//     for (int counter=1;counter<NumEntries;counter++){
//       CycleTable(counter).C=CycleTable(counter-1).C+CycleTable(counter).P;
//     }
//   }
//   else if (OnlyEven) {
//     int currSpot=0;
//     for (int counter=0;counter<NumEntries;counter++){
//       if (CycleTable(counter).Length % 2 ==0){
// 	CycleTable(currSpot)=CycleTable(counter);
// 	currSpot++;
//       }
//     }
//     NumEntries=currSpot;
//     CycleTable(0).C=CycleTable(0).P;
//     for (int counter=1;counter<NumEntries;counter++){
//       CycleTable(counter).C=CycleTable(counter-1).C+CycleTable(counter).P;
//     }
//   }

  Norm = CycleTable(NumEntries-1).C;
  NormInv = 1.0/Norm;
  //  cerr<<"The number of my table is "<<NumEntries<<endl;
}
