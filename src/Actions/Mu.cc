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

#include "../PathDataClass.h"
#include "Mu.h"

///This has to be called after pathdata knows how many
///particles it has
void MuClass::Read(IOSectionClass& in)
{
  in.ReadVar("Mu",Mu);
  //  Mu=-0.72/10+0.9633;
  Mu=1.0;
}

MuClass::MuClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
}

bool 
MuClass::PadWorm()
{
#if 1==2
  int headSlice;
  int headPtcl;
  int tailSlice;
  int tailPtcl;
  PathData.FindHead(headSlice,headPtcl);
  PathData.FindTail(tailSlice,tailPtcl);


  int nextPtcl=tailPtcl;
  int nextSlice=tailSlice;
  int numEmpty=0;
  do{
    if (nextSlice==PathData.Join)
      nextPtcl=PathData.Path.Permutation(nextPtcl);
    nextSlice=(nextSlice+1) % PathData.Path.NumTimeSlices();
    if (nextSlice!=0) //We should only count iether the beginning or the end slice but not both of them
      numEmpty++;
  }
  while (PathData.Path.ParticleExist(nextSlice,nextPtcl)==0.0);
  numEmpty--;
  cerr<<"Num full is "<<PathData.Path.NumTimeSlices()*PathData.Path.NumParticles()-numEmpty-50<<endl;
  if (numEmpty<PathData.Path.NumTimeSlices()+1){
    return true;
  }
  else 
    return false;
//   if (numEmpty<PathData.Path.NumParticles()){
//     int oldNumParticles=PathData.Path.NumParticles();
//     int newNumParticles=oldNumParticles+1;
//     PathData.Path.resizeAndPreserve(PathData.Path.NumTimeSlices(),
// 				    newNumParticles);
//     PathData.Path.ParticleExist.resizeAndPreserve(PathData.Path.NumTimeSlices(),
// 						  newNumParticles);
//     PathData.Path.Permutation.resizeAndPreserve(newNumParticles);
//   }
#endif
}

double 
MuClass::SingleAction (int slice1, int slice2,
			       const Array<int,1> &changedParticles,
			       int level)
{
  //  cerr<<"In and out of mu"<<endl;
  return 0.0;
#if 1==2
  int headPtcl;
  int tailPtcl;
  int headSlice;
  int tailSlice;
  double lambda;
  PathData.FindHead(headSlice,headPtcl);
  PathData.FindTail(tailSlice,tailPtcl);
  //  cerr<<headSlice<<" "<<headPtcl<<" "<<tailSlice<<" "<<tailPtcl<<endl;
  //  PathData.Path.PrintRealSlices();
//   double toReturn=0.0;
//   if (PathData.Path.NowOpen){
//     dVec disp=PathData.Path(headSlice,headPtcl)-PathData.Path(tailSlice,tailPtcl);
//     PathData.Path.PutInBox(disp);
//     double dist2=dot(disp,disp);
//     lambda=
//       PathData.Path.Species(PathData.Path.ParticleSpeciesNum(headPtcl)).lambda;
//     toReturn+=dist2/(4.0*lambda*PathData.Path.tau);
//   }
//   else {
//     toReturn+=0.0;
//   }
    

  double totalMu=0;
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
  //  cerr << "I'm returning kinetic action " << TotalK << endl;

  int numSlices=0;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
      totalMu+=PathData.Path.ParticleExist(slice,ptcl);
      if (PathData.Path.ParticleExist(slice,ptcl)==1.0)
	numSlices++;
    }
  totalMu*=Mu/(double)PathData.Path.NumTimeSlices();
    int wormSize,numEmpty;
    //  int tailSlice,tailPtcl;
    //  int headSlice,headPtcl;
//    PathData.WormInfo(headSlice, headPtcl,
//  		    tailSlice, tailPtcl,
//  		    numEmpty, wormSize);
  

    int numFull=0;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
      if (PathData.Path.IsFull(ptcl)){
	numFull++;
      }
  
  if (!PathData.Path.NowOpen){
    //    double numParticles=(double)wormSize/(double)PathData.Path.NumTimeSlices();
    cerr<<"NuMParticles: "<<numFull<<" "<<wormSize<<" "<<PathData.Path.NumTimeSlices()<<endl;
    if (numFull!=0)
      totalMu=totalMu-log(5.0)-log((double)PathData.Path.NumTimeSlices()-1)-1.0*log((double)numFull);
  }
    
  //  totalMu+=1.5*(double)numSlices*log(4*M_PI*lambda*PathData.Path.tau);
//   if (PadWorm())
//     return -100000;
//  else
  //return totalMu;
#endif
}


string
MuClass::GetName()
{
  return "Mu";
}
