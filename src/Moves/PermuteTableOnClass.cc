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

// #include "PermuteTableClass.cc"
#include "PermuteTableClass.h"
#include "PermuteTableOnClass.h"


double PermuteTableOnClass::AttemptPermutation()
{
  //Get a random number number from the local processor stream.
//   Array<int,1> TestArray(NumEntries);
//   TestArray=0;
//   int index;
//   for (int counter=0;counter<1000000;counter++){
//     double xi=PathData.Path.Random.Local(); 
//     index=FindEntry(xi);
//     TestArray(index)++;
//   }
//   //  PrintTable();
//   //  cerr<<"Printing stuff: "<<endl;
//   for (int counter=0;counter<TestArray.size();counter++){
//     cerr<<TestArray(counter)<<endl;
//   }
  double xi=PathData.Path.Random.Local(); 
  int index=FindEntry(xi);
  CurrentCycle = CycleTable(index);
  //  cerr<<"Index is "<<xi<<" "<<index<<" "<<endl;
  //  for (int j=0; j<CurrentCycle.Length; j++)
  //    cerr << CurrentCycle.CycleRep[j] << " ";
  //  cerr<<"Done"<<endl;
//   if (CurrentCycle.CycleRep[0]<0 || CurrentCycle.CycleRep[0]>10000){
//     PrintTable();
//   }
  //  cerr<<"After find entry"<<endl;
  // Now, apply the permutation to the Path
  int firstPtcl=PathData.Species(SpeciesNum).FirstPtcl;
  CurrentCycle.Apply(PathData.Path,firstPtcl,Slice2);
  
//   int diff = Slice2-Slice1;
//   int numLevels = 0;
//   while (diff != 1) {
//     diff >>= 1;
//     numLevels++;
//   }
  ///Remember to update distance table after permutation


  return (CurrentCycle.P * NormInv);
}

// double PermuteTableOnClass::Sample (int &slice1, int &slice2,
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

double PermuteTableOnClass::CalcReverseProb(const PermuteTableOnClass &forwardTable)
{
  
  // Reverse probabily for a single ptcl move is the same as that for the
  // forward move.
  if (forwardTable.CurrentCycle.Length == 1)
    return (forwardTable.CurrentCycle.P*forwardTable.NormInv);
  //We reconstruct things from scratch to make sure we do the right
  //thing. We can try to incrementally make this faster.
  ConstructCycleTable(forwardTable.SpeciesNum,
		      forwardTable.Slice1,forwardTable.Slice2,
		      forwardTable.ExcludeParticle); 
  //  cerr << "Forward Table:\n";
  //  forwardTable.PrintTable();
  //  cerr << "Reverse Table:\n";
  //  PrintTable();
  ///Correct for gamma factors here


//   int i=0;
//   while (!(CycleTable(i)==forwardTable.CurrentCycle))
//     i++;
//   if (fabs(CycleTable(i).P-
// 	   (1.0/(forwardTable.CurrentCycle.P)))>1e-12)
//     cerr <<"BAD! BAD! BAD! "
// 	 << CycleTable(i).P<<" "
// 	 <<1.0/(forwardTable.CurrentCycle.P);
//   }
  
  return (1.0/(forwardTable.CurrentCycle.P)*NormInv);
	 
	   
}
    
Array<int,1> PermuteTableOnClass::CurrentParticles()
{
  Array<int,1> tempPtcl(CurrentCycle.Length);
  int firstPtcl=PathData.Species(SpeciesNum).FirstPtcl;
  for (int i=0;i<CurrentCycle.Length;i++){
    tempPtcl(i)=CurrentCycle.CycleRep(i)+firstPtcl;
  }
  return tempPtcl;
}

double PermuteTableOnClass::calcij(int i, int j)
{
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Path.tau * (double) (Slice2-Slice1);
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  int N = lastPtcl-firstPtcl+1;
  dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
  double dist_ii = dot (disp_ii,disp_ii);
  PathData.Path.PutInBox(disp_ii);
  dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
  PathData.Path.PutInBox(disp_ij);
  double dist_ij = dot (disp_ij,disp_ij);
  return exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
}

void PermuteTableOnClass::ConstructHTable()
{
  PathClass &Path=PathData.Path;
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Path.tau * (double) (Slice2-Slice1);
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  int N = lastPtcl-firstPtcl+1;
  for (int i=0;i<N;i++){
    dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
    PathData.Path.PutInBox(disp_ii);
    double dist_ii = dot (disp_ii,disp_ii);
    int xBox,yBox,zBox;
    ////BUG!    PathData.Path.Cell.FindBox(PathData.Path(Slice1,i+firstPtcl),xBox,yBox,zBox);
    int rxbox,rybox,rzbox;
    for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
      rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
      rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
      rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
      list<int> &ptclList=PathData.Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(Slice2);
      for (list<int>::iterator iter=ptclList.begin();iter!=ptclList.end();iter++) {
	int j=*iter;
	dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
	PathData.Path.PutInBox(disp_ij);
	double dist_ij = dot (disp_ij,disp_ij);
	ParticleProbStruct particleProb;
	particleProb.ptcl=j;
	particleProb.prob=exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
	if (particleProb.prob>epsilon)
	  HTable(i).push_back(particleProb);
      }
    }
  }
}

 
///// This constructs the Htable for the reverse move from the Htable
///// for the forward move.  i.e. Constructs Hrev(i,j) = H(P(i),P(j)).
// void PermuteTableOnClass::PermuteHTable(const CycleClass &myPerm,
// 				      const PermuteTableOnClass &forwardTable)
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
// void CycleClass::Apply(PathClass &path, int firstPtcl, int slice)
// {
//   SetMode(NEWMODE);;


void PermuteTableOnClass::ConstructCycleTable(int speciesNum,int slice1,int slice2)
{
  ConstructCycleTable(speciesNum,slice1,slice2,-1);
}


void PermuteTableOnClass::ConstructCycleTable(int speciesNum,
					    int slice1, int slice2,
					    int particleExclude)
{
  Slice1=slice1;
  Slice2=slice2;
  SpeciesNum=speciesNum;
  ExcludeParticle=particleExclude;
  ConstructHTable();
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  CycleTable.resize(TableSize);
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
      tempPerm.C=tempPerm.P;
      if (NumEntries!=0){
	tempPerm.C+=CycleTable(NumEntries-1).C;
      }
      AddEntry(tempPerm);
      //      for (int j=i+1; j<N; j++) {// 2 and higher cycles
      for (list<ParticleProbStruct>::iterator jIter=HTable(i).begin();
	   jIter!=HTable(i).end();jIter++) {
	int j=(*jIter).ptcl;
	if (j!=particleExclude && i!=particleExclude && j>i){
	  //2 cycle permutations
	  tempPerm.CycleRep[1]=j;
	  hprod=(*jIter).prob;
	  tempPerm.Length=2;
	  tempPerm.P=Gamma(1)*hprod*calcij(j,i);
	  tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;     
	  if (tempPerm.P > epsilon)
	    AddEntry(tempPerm);
	  if (hprod * Gamma(2)*Gamma(3)>epsilon) {
	    //	    for (int k=i+1;k<N;k++){//3 and higher cycles
	    for (list<ParticleProbStruct>::iterator kIter=HTable(j).begin();
		 kIter!=HTable(j).end();kIter++) {
	      int k=(*kIter).ptcl;
	      if (k!=particleExclude && k!=j && k>i){
		//3 cycle permutations
		tempPerm.CycleRep[2]=k;
		hprod2=hprod*(*kIter).prob; //HTable(j,k);
		tempPerm.Length=3;
		tempPerm.P=Gamma(2)*hprod2*calcij(k,i); //HTable(k,i);
		tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;
		if (tempPerm.P > epsilon)
		  AddEntry(tempPerm);
		if (hprod2 *Gamma(3) > epsilon) {
		  for (list<ParticleProbStruct>::iterator lIter=HTable(k).begin();
		       lIter!=HTable(k).end();lIter++) {
		    int l=(*lIter).ptcl;
		    if (l!=particleExclude && (l!=j) && (l!=k) && l>i){
		      hprod3=hprod2*(*lIter).prob; //HTable(k,l);
		      tempPerm.CycleRep[3]=l;
		      tempPerm.Length=4;
		      tempPerm.P=Gamma(3)*hprod3*calcij(l,i); //HTable(l,i);
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
  if (OnlyOdd || PathData.Path.Species(SpeciesNum).GetParticleType()==FERMION){
    int currSpot=0;
    for (int counter=0;counter<NumEntries;counter++){
      if (CycleTable(counter).Length % 2 ==1){
	CycleTable(currSpot)=CycleTable(counter);
	currSpot++;
      }
    }
    NumEntries=currSpot;
    CycleTable(0).C=CycleTable(0).P;
    for (int counter=1;counter<NumEntries;counter++){
      CycleTable(counter).C=CycleTable(counter-1).C+CycleTable(counter).P;
    }
  }
  else if (OnlyEven) {
    int currSpot=0;
    for (int counter=0;counter<NumEntries;counter++){
      if (CycleTable(counter).Length % 2 ==0){
	CycleTable(currSpot)=CycleTable(counter);
	currSpot++;
      }
    }
    NumEntries=currSpot;
    CycleTable(0).C=CycleTable(0).P;
    for (int counter=1;counter<NumEntries;counter++){
      CycleTable(counter).C=CycleTable(counter-1).C+CycleTable(counter).P;
    }
  }

  Norm = CycleTable(NumEntries-1).C;
  NormInv = 1.0/Norm;
  //  cerr<<"The number of my table is "<<NumEntries<<endl;
}




// ///the Gamma's must all be greater then 1 or this won't do the correct
// ///thing.  Allows you to exclude a particle when producing the table
// ///(say for example the open particle).
// void PermuteTableOnClass::ConstructCycleTable(int speciesNum,
// 					    int slice1, int slice2,
// 					    int particleExclude)
// {
//   Slice1=slice1;
//   Slice2=slice2;
//   SpeciesNum=speciesNum;
//   ExcludeParticle=particleExclude;
//   ConstructHTable();
//   int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
//   int lastPtcl= PathData.Species(SpeciesNum).LastPtcl;
//   CycleTable.resize(TableSize);
//   NumEntries = 0;
//   int N=HTable.extent(0);
//   int xBox,yBox,zBox;
//   int xEffect=PathData.Path.Cell.Xeffect;
//   int yEffect=PathData.Path.Cell.Yeffect;
//   int zEffect=PathData.Path.Cell.Zeffect;
//   double lambda = PathData.Species(SpeciesNum).lambda;
//   double beta = PathData.Path.tau * (double) (Slice2-Slice1);
//   double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
//   for (int ptcl1=firstPtcl;ptcl1<lastPtcl;ptcl1++){
//     CycleClass tempPerm;
//     tempPerm.Length=1;
//     tempPerm.CycleRep[0]=ptcl1;
//     tempPerm.P=Gamma(0);
//     tempPerm.C=tempPerm.P;
//     if (NumEntries!=0){
//       tempPerm.C+=CycleTable(NumEntries-1).C;
//     }
//     AddEntry(tempPerm);
//     PathData.Path.Cell.FindBox(PathData.Path(Slice1,ptcl1),xBox,yBox,zBox);
//     for (int xCell=xBox-xEffect;xCell<=xBox+xEffect;xCell++){
//       for (int yCell=yBox-yEffect;yCell<=yBox+yEffect;yCell++){
// 	for (int zCell=zBox-zEffect;zCell<=zBox+zEffect;zCell++){
// 	  int rxbox,rybox,rzbox;
// 	  rxbox=(xCell+PathData.Path.Cell.GridsArray.extent(0)) % PathData.Path.Cell.GridsArray.extent(0);
// 	  rybox=(yCell+PathData.Path.Cell.GridsArray.extent(1)) % PathData.Path.Cell.GridsArray.extent(1);
// 	  rzbox=(zCell+PathData.Path.Cell.GridsArray.extent(2)) % PathData.Path.Cell.GridsArray.extent(2);
// 	  //	    cerr<<rxbox<<" "<<rybox<<" "<<rzbox<<endl;
// 	  list<int> &ptclList=PathData.Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(Slice2);
// 	  for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
// 	    int ptcl2=*i;
// 	    if (ptcl2!=ptcl1){
// 	      tempPerm.Length=2;
// 	      tempPerm.CycleRep[1]=ptcl2;
	      
// 	      dVec disp_ii = PathData(Slice2,ptcl1)-PathData(Slice1,ptcl1);
// 	      PathData.Path.PutInBox(disp_ii);
// 	      double dist_ii = dot (disp_ii,disp_ii);
// 	      dVec disp_ij = PathData(Slice2,ptcl2)-PathData(Slice1,ptcl1);
// 	      PathData.Path.PutInBox(disp_ij);
// 	      double dist_ij = dot (disp_ij,disp_ij);
// 	      double dist1 = exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
	      
	      
// 	      disp_ii = PathData(Slice2,ptcl2)-PathData(Slice1,ptcl2);
// 	      PathData.Path.PutInBox(disp_ii);
// 	      dist_ii = dot (disp_ii,disp_ii);
// 	      disp_ij = PathData(Slice2,ptcl1)-PathData(Slice1,ptcl2);
// 	      PathData.Path.PutInBox(disp_ij);
// 	      dist_ij = dot (disp_ij,disp_ij);
// 	      double dist2 = exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
	      
	      
// 	      tempPerm.P=Gamma(1)*dist1*dist2;
// 	      tempPerm.C=tempPerm.P+CycleTable(NumEntries-1).C;     
// 	      if (tempPerm.P > epsilon)
// 		AddEntry(tempPerm);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
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

//   Norm = CycleTable(NumEntries-1).C;
//   NormInv = 1.0/Norm;
//   //  cerr<<"The number of my table is "<<NumEntries<<endl;
// }
