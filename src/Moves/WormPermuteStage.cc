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
#include "WormPermuteStage.h"

double WormPermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}

void WormPermuteStageClass::Accept()
{
  NumAccepted(0)++;
  NumAttempted(0)++;
}

void WormPermuteStageClass::Reject()
{
  NumAttempted(0)++;
}


void WormPermuteStageClass::WriteRatio()
{
//   Array<int,1> numAttemptTotal(4);
//   Array<int,1> numAcceptTotal(4);
//   Array<double,1> ratioTotal(4);
//   int totalAttempts=0;
//   cerr<<"A"<<endl;
//   PathData.Path.Communicator.Sum(NumAttempted,numAttemptTotal);
//   cerr<<"A2"<<endl;
//   PathData.Path.Communicator.Sum(NumAccepted,numAcceptTotal);
//   for (int len=0;len<4;len++){
//     totalAttempts=totalAttempts+numAcceptTotal(len);
//     if (numAttemptTotal(len)!=0)

//       ratioTotal(len)=(double)numAcceptTotal(len)/((double)numAttemptTotal(len));
//     else
//       ratioTotal(len)=0.0;
//   }
//   cerr<<"B"<<endl;
//       ///divides by 2because accept gets called twice in acepting stages
//   for (int i=0;i<numAttemptTotal.size();i++)
//     numAttemptTotal(i)=(int)(numAttemptTotal(i)/2);
//   cerr<<"C"<<endl;
 
//   if (totalAttempts!=0){
//     AcceptanceRatioVar.Write(ratioTotal);
//     AcceptanceTotalVar.Write(numAttemptTotal);
//     NumAttempted=0;
//     NumAccepted=0;
//   }
//   cerr<<"out of my permute stage ratio"<<endl;

}



void WormPermuteStageClass::Read (IOSectionClass &in)
{
  OnlyX=false;
  in.ReadVar("OnlyX",OnlyX);
  LocalStageClass::Read(in);
  Initialize();
}

//Decides which particles to attempt to permute onto. If you have set
//OnlyX then it  choose to permute onto particles that are started
//with the x value within 0.01 of the origin. Otherwise it attempts to
//permute onto all particles.  
//Must be called after the path has been initalized.
void WormPermuteStageClass::Initialize()
{
  cerr<<"I've started initializing"<<endl;
  if (OnlyX){
    cerr<<"The number of particles is "<<PathData.Path.NumParticles()
	<<PermuteOntoList.size()<<endl;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      if (abs(PathData.Path(0,ptcl)[1]-0.0)<0.01 && 
	  abs(PathData.Path(0,ptcl)[2]-0.0)<0.01){
	PermuteOntoList.push_back(ptcl);
      }
    }
  }
  else {
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      PermuteOntoList.push_back(ptcl);
    }
  }
    
  cerr<<"The size of my list is "<<PermuteOntoList.size()<<endl;
}

int WormPermuteStageClass::CountPtclInOpenLoop()
{
  int currPtcl=PathData.Path.OpenPtcl;
  int totalCount=1;
  while (PathData.Path.Permutation(currPtcl)!=PathData.Path.OpenPtcl){
    totalCount++;
    currPtcl=PathData.Path.Permutation(currPtcl);
  }
  return totalCount;
}

bool WormPermuteStageClass::PtclInOpenLoop(int checkPtcl)
{
  int currPtcl=PathData.Path.OpenPtcl;
  if (currPtcl==checkPtcl)
    return true;
  while (PathData.Path.Permutation(currPtcl)!=PathData.Path.OpenPtcl){
    currPtcl=PathData.Path.Permutation(currPtcl);
    if (currPtcl==checkPtcl)
      return true;
  }
  return false;
  

}

void WormPermuteStageClass::InitBlock(int &slice1,int &slice2)
{


  ///Move the path so that it's setup with the tail at the last time slice
  while (PathData.Path.OpenLink!=PathData.Path.NumTimeSlices()-1){
    PathData.MoveJoin(0);
    PathData.ShiftData(1);
    PathData.Join=1;
  }
  
  PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
  slice2=PathData.Path.NumTimeSlices()-1;
  slice1=0;
  //  slice2=slice1+(1<<NumLevels);
  
  if (PathData.Path.OrderN){
    for (int slice=slice1;slice<=slice2;slice++)
      PathData.Path.Cell.BinParticles(slice);
  }
  
}

bool WormPermuteStageClass::Attempt (int &slice1, int &slice2, 
				      Array<int,1> &activeParticles, double &prevActionChange)
{

//   //sortof hackish here
//   slice1=Slice1;
//   slice2=Slice2;
//  cerr<<"Inside attempt"<<endl;
  if (activeParticles(0)!=-1)
    return true;
  int numOfPtclInOpenLoop=CountPtclInOpenLoop();
//   if (numOfPtclInOpenLoop>=3){
//     cerr<<"end info"<<endl;
//     cerr<<"Info: "<<endl;
//     cerr<<PathData.Path.OpenLink<<endl;
//     cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.OpenPtcl)<<endl;
//     cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())<<endl;
//     cerr<<PathData.Path.OpenPtcl;
//     cerr<<numOfPtclInOpenLoop<<endl;
//     cerr<<slice1<<" "<<slice2<<endl;
//     for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
//       cerr<<PathData.Path.Permutation(ptcl)<<" ";
//     cerr<<endl;
//     assert(1==2);
//   }
  int permuteWith=PathData.Path.OpenPtcl;
  while (permuteWith==PathData.Path.OpenPtcl)
    permuteWith=
      PermuteOntoList[PathData.Path.Random.LocalInt(PermuteOntoList.size())];
  activeParticles.resize(2);
  
  
  if (!PtclInOpenLoop(permuteWith)){

      activeParticles(0)=permuteWith;
      activeParticles(1)=PathData.Path.OpenPtcl;
      dVec currTail=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.NumParticles());
      dVec currHead=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.OpenPtcl);
      int lastSlice=PathData.NumTimeSlices()-1;
      int tailParticle=PathData.Path.NumParticles();
      int openPtcl=PathData.Path.OpenPtcl;
      PathData.Path(lastSlice,tailParticle)=PathData.Path(lastSlice,permuteWith);
      PathData.Path(lastSlice,openPtcl)=currTail;
      PathData.Path(slice1,permuteWith)=currTail;
      PathData.Path(lastSlice,permuteWith)=currHead;

      //Swap permutation
      int oldPerm=PathData.Path.Permutation(PathData.Path.OpenPtcl);
      PathData.Path.Permutation(PathData.Path.OpenPtcl)=PathData.Path.Permutation(permuteWith);
      PathData.Path.Permutation(permuteWith)=oldPerm;
      //end swap
      PathData.Path.OpenPtcl=permuteWith;
      prevActionChange=0;
      int numOfPtclInOpenLoop=CountPtclInOpenLoop();
      //      cerr<<"The number of particles in the open loop are A "<<numOfPtclInOpenLoop<<endl;
      //      cerr<<"Leaving worm stage at this spot"<<endl
//       if (numOfPtclInOpenLoop>=2){
// 	cerr<<"Info: "<<endl;
// 	cerr<<PathData.Path.OpenLink<<endl;
// 	cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.OpenPtcl)<<endl;
// 	cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())<<endl;
// 	cerr<<PathData.Path.OpenPtcl;
// 	cerr<<permuteWith<<endl;
// 	cerr<<numOfPtclInOpenLoop<<endl;
// 	cerr<<slice1<<" "<<slice2<<endl;
//       }
      return true;
    }
    else{
      int permutedOntoOpen=-1;
      for (int ptcl=0;ptcl<PathData.NumParticles();ptcl++)
	if (PathData.Path.Permutation(ptcl)==PathData.Path.OpenPtcl)
	  permutedOntoOpen=ptcl;
      activeParticles(0)=permutedOntoOpen;
      activeParticles(1)=PathData.Path.OpenPtcl;

      dVec currTail=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.NumParticles());
      dVec currHead=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.OpenPtcl);
      int lastSlice=PathData.NumTimeSlices()-1;
      int tailParticle=PathData.Path.NumParticles();
      int openPtcl=PathData.Path.OpenPtcl;
      PathData.Path(lastSlice,openPtcl)=currTail;
      PathData.Path(slice1,openPtcl)=currTail;
      PathData.Path(lastSlice,tailParticle)=
	PathData.Path(lastSlice,permutedOntoOpen);
      PathData.Path(lastSlice,permutedOntoOpen)=currHead;

      //Swap permutation
      int oldPerm=PathData.Path.Permutation(PathData.Path.OpenPtcl);
      PathData.Path.Permutation(PathData.Path.OpenPtcl)=PathData.Path.Permutation(permutedOntoOpen);
      PathData.Path.Permutation(permutedOntoOpen)=oldPerm;
      //end swap
      
      PathData.Path.OpenPtcl=permutedOntoOpen;
      prevActionChange=0;

      int numOfPtclInOpenLoop=CountPtclInOpenLoop();
//       if (numOfPtclInOpenLoop>=2){
// 	cerr<<"Info Other: "<<endl;
// 	cerr<<PathData.Path.OpenLink<<endl;
// 	cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.OpenPtcl)<<endl;
// 	cerr<<PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())<<endl;
// 	cerr<<PathData.Path.OpenPtcl;
// 	cerr<<permutedOntoOpen<<endl;
// 	cerr<<numOfPtclInOpenLoop<<endl;
//       }
      //      cerr<<"The number of particles in the open loop are "<<numOfPtclInOpenLoop<<endl;
      //      cerr<<"Leaving worm stage at this other spot"<<endl;
      return (PathData.Path.Random.Local()<(double)(1.0/(double)(numOfPtclInOpenLoop)));
	//       if  (numOfPtclInOpenLoop==1)
// 	return true;
//       else if (numOfPtclInOpenLoop==2)
// 	return (PathData.Path.Random.Local()<1.0/2.0);
//       else if (numOfPtclInOpenLoop==3)
// 	return (PathData.Path.Random.Local()<1.0/3.0);
//       else if (numOfPtclInOpenLoop==4)
// 	return (PathData.Path.Random.Local()<1.0/4.0);
      //}
    }

}



// bool WormPermuteStageClass::Attempt (int &slice1, int &slice2, 
// 				      Array<int,1> &activeParticles, double &prevActionChange)
// {
//   if (activeParticles(0)==-1){

//     //    while (PathData.Path.OpenLink-(1<<NumLevels)>=1 && PathData.Path.OpenLink<=PathData.Path.NumTimeSlices()-1)){
//     while (PathData.Path.OpenLink!=PathData.Path.NumTimeSlices()-1){
//       PathData.MoveJoin(0);
//       PathData.ShiftData(1);
//       PathData.Join=1;
//     }

//     PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
//     slice1=0;
//     slice2=PathData.Path.NumTimeSlices()-1;

//     int permuteWith=PathData.Path.OpenPtcl;
//     while (permuteWith==PathData.Path.OpenPtcl)
//       permuteWith=PathData.Path.Random.LocalInt(9);
//     activeParticles.resize(2);

// //     cerr<<"My current loops is ";
// //     int ptclLoop=PathData.Path.OpenPtcl;
// //     cerr<<ptclLoop<<" ";
// //     while (PathData.Path.Permutation(ptclLoop)!=PathData.Path.OpenPtcl){
// //       ptclLoop=PathData.Path.Permutation(ptclLoop);
// //       cerr<<ptclLoop<<" ";
// //     }
// //     cerr<<endl;
// //     cerr<<"And trying to permute with "<<permuteWith<<endl;

//     if (PathData.Path.Permutation(permuteWith)==permuteWith){

//     activeParticles(0)=permuteWith;
//     activeParticles(1)=PathData.Path.OpenPtcl;

//     dVec currTail=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.NumParticles());

//     dVec currHead=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.OpenPtcl);
    
//     int lastSlice=PathData.NumTimeSlices()-1;
//     int tailParticle=PathData.Path.NumParticles();
//     int openPtcl=PathData.Path.OpenPtcl;

//     PathData.Path(lastSlice,tailParticle)=PathData.Path(lastSlice,permuteWith);
//     PathData.Path(lastSlice,openPtcl)=currTail;
//     PathData.Path(slice1,permuteWith)=currTail;
    
//     PathData.Path(lastSlice,permuteWith)=currHead;


//     //Swap permutation
//     int oldPerm=PathData.Path.Permutation(PathData.Path.OpenPtcl);
//     PathData.Path.Permutation(PathData.Path.OpenPtcl)=PathData.Path.Permutation(permuteWith);
//     PathData.Path.Permutation(permuteWith)=oldPerm;
//     //end swap


//     PathData.Path.OpenPtcl=permuteWith;
//     prevActionChange=0;
//     //    if (PathData.Path.Permutation(0)!=0 && 
//     //	PathData.Path.Permutation(1)!=1 &&
//     //	PathData.Path.Permutation(2)!=2)
//     //      return (PathData.Path.Random.Local()>0.5);
//     //    else
//     return true;
//     }
//     else{
//       dVec currTail=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.NumParticles());

//       dVec currHead=PathData.Path(PathData.Path.NumTimeSlices()-1,PathData.Path.OpenPtcl);
    
//       int lastSlice=PathData.NumTimeSlices()-1;
//       int tailParticle=PathData.Path.NumParticles();
//       int openPtcl=PathData.Path.OpenPtcl;
//       PathData.Path(lastSlice,openPtcl)=currTail;
//       PathData.Path(slice1,openPtcl)=currTail;
//       int permutedOntoOpen=-1;
//       for (int ptcl=0;ptcl<PathData.NumParticles();ptcl++)
// 	if (PathData.Path.Permutation(ptcl)==PathData.Path.OpenPtcl)
// 	  permutedOntoOpen=ptcl;
//       PathData.Path(lastSlice,tailParticle)=
// 	PathData.Path(lastSlice,permutedOntoOpen);
//       PathData.Path(lastSlice,permutedOntoOpen)=currHead;

//     //Swap permutation
//     int oldPerm=PathData.Path.Permutation(PathData.Path.OpenPtcl);
//     PathData.Path.Permutation(PathData.Path.OpenPtcl)=PathData.Path.Permutation(permutedOntoOpen);
//     PathData.Path.Permutation(permutedOntoOpen)=oldPerm;
//     //end swap
      
//     activeParticles(0)=permutedOntoOpen;
//     activeParticles(1)=PathData.Path.OpenPtcl;

//     PathData.Path.OpenPtcl=permutedOntoOpen;
//     prevActionChange=0;
//     //    if (PathData.Path.Permutation(PathData.Path.OpenPtcl)==PathData.Path.OpenPtcl)
//     //      return (PathData.Path.Random.Local()>0.5);
//       //    else 
//     if (PathData.Path.Permutation(PathData.Path.OpenPtcl)==PathData.Path.OpenPtcl)
//       return true;
//     else
//       return (PathData.Path.Random.Local()>0.5);

//     }
//   }
//   else{
  
  
//     return true;    
//   }
// }
