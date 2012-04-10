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

#include "RandomPermClass.h"

void  RandomPermClass::GeneratePerm()
{
  firstIndex i;
  PermArray=i;
  int numPtcl=PermArray.size();
  for (int mcStep=0;mcStep<numPtcl*numPtcl;mcStep++){
    int l1=(int)(PathData.Path.Random.Local()*numPtcl);
    int l2=(int)(PathData.Path.Random.Local()*numPtcl);
    //    cerr<<(ActionPairs(l2,PermArray(l1))*ActionPairs(l1,PermArray(l2)))/(ActionPairs(l1,PermArray(l1))*ActionPairs(l2,PermArray(l2)))<<endl;
    if ((ActionPairs(l2,PermArray(l1))*ActionPairs(l1,PermArray(l2)))/(ActionPairs(l1,PermArray(l1))*ActionPairs(l2,PermArray(l2)))>
	PathData.Path.Random.Local()){
      int temp=PermArray(l1);
      PermArray(l1)=PermArray(l2);
      PermArray(l2)=temp;
    }
  }
  /////  for (int counter=0;counter<PermArray.size();counter++){
  ////    cerr<<PermArray(counter)<<" ";
  ////}
/////  cerr<<endl;
  return;
}
void RandomPermClass::Apply()
{
  int totalPtcl=0;
  double probRatio=0;
  int firstPtcl=PathData.Path.Species(SpeciesNum).FirstPtcl;
  for (int counter=0;counter<PermArray.size();counter++){
    if (PermArray(counter) != counter){
      totalPtcl++;
      SetMode(OLDMODE);
      dVec newPlace=PathData.Path(Slice2,PermArray(counter));
      SetMode(NEWMODE);
      PathData.Path.SetPos(Slice2,counter,newPlace);
    } 
  }
  int currSpot=0;
  CurrentParticles.resize(totalPtcl);
  for (int counter=0;counter<PermArray.size();counter++){
    if (PermArray(counter)!=counter){
      CurrentParticles(currSpot)=counter;
      currSpot=currSpot+1;
      probRatio += log(ActionPairs(counter,PermArray(counter)))-log(ActionPairs(counter,counter));
    }
  }
  
}
void RandomPermClass::GetActionPairs()
{
  double beta = PathData.Path.tau * (double) (Slice2-Slice1);
  double lambda = PathData.Species(SpeciesNum).lambda;
  //  cerr<<"Beta is "<<beta<<endl;

  //  cerr<<"Lambda is "<<lambda<<endl;
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  //  cerr<<"Four lambda beta is "<<fourLambdaBetaInv<<endl;
  int firstPtcl=PathData.Path.Species(SpeciesNum).FirstPtcl;
  int lastPtcl=PathData.Path.Species(SpeciesNum).LastPtcl;
  dVec disp_ij;
  double dist_ij;


  for (int ptcl1=firstPtcl;ptcl1<=lastPtcl;ptcl1++){
    for (int ptcl2=firstPtcl;ptcl2<=lastPtcl;ptcl2++){
      disp_ij=PathData.Path(Slice2,ptcl2)-PathData.Path(Slice1,ptcl1); //things break if landing on open
      PathData.Path.PutInBox(disp_ij);
      dist_ij=dot(disp_ij,disp_ij);
      ActionPairs(ptcl1,ptcl2)=exp(-dist_ij*fourLambdaBetaInv);
    }
  }
  
  
}


void RandomPermClass::Read(IOSectionClass &inSection)
{

  
  //  Array<double,1> tempGamma;
  //  assert(inSection.ReadVar("Gamma",tempGamma));
//  assert(tempGamma.size() == 4);
//  for (int counter=0;counter<tempGamma.size();counter++){
//    Gamma(counter)=tempGamma(counter);
// }
//  assert(inSection.ReadVar("epsilon",epsilon));
  //  assert(inSection.ReadVar("SpeciesNum",SpeciesNum));
 
  
}

RandomPermClass::RandomPermClass(PathDataClass &myPathData) : PathData(myPathData)
{
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  int N = lastPtcl-firstPtcl+1;
  ActionPairs.resize(N,N);
  PermArray.resize(N);
}
