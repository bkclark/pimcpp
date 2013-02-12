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

#ifndef RANDOM_PERM_H
#define RANDOM_PERM_H

#include "MoveBase.h"
#include "RandomPermClass.h"

class PermMove : public ParticleMoveClass
{
//  public:
//   RandomPermClass myPerm;
//   BisectionClass Bisection;
//   int NumLevels;
//   void MakeMove()
//   {
//     // First, decide on the chunk of slices we're working on
//     int slice1;
//     int slice2;
//     int numSlices=PathData.NumTimeSlices();
    
//     if (!PathData.Path.OpenPaths){
//       double xi=PathData.Path.Random.Local();
//       slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
//       slice1=slice2-(1<<NumLevels);
//     }
//     else {
//       // First, decide on the chunk of slices we're working on
//       do{
// 	double xi=PathData.Path.Random.Local();
// 	slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
// 	slice1=slice2-(1<<NumLevels);
// 	//      if (slice2>=PathData.Path.NumTimeSlices()-1){
// 	//	cerr<<"ERROR!"<<slice2<<endl;
// 	//	slice2--;
// 	//      }
//       } while ((slice1<=(int)(PathData.Path.OpenLink) && (int)(PathData.Path.OpenLink)<=slice2));
//     }
//     //  if (slice2>=PathData.Path.NumTimeSlices()-1){
//     //    cerr<<"ERROR! ERROR! ERROR!";
//     //  }
//     //  cerr<<slice1<<" "<<slice2<<" "<<PathData.Path.NumTimeSlices()<<endl;
//     PathData.MoveJoin(slice2);
//     myPerm.SpeciesNum=0;
//     myPerm.Slice1=slice1;
//     myPerm.Slice2=slice2;
//     myPerm.GetActionPairs();

//     for (int counter=0;counter<100;counter++){
//       myPerm.GeneratePerm();
//       //    firstIndex i;
//       //    Array<int,1> currentParticles(PathData.Path.NumParticles());
//       //    currentParticles=i;
//       myPerm.Apply();
//       bool acceptBisect;
//       if (myPerm.CurrentParticles.size()!=0){
// 	acceptBisect = 
// 	  Bisection.Bisect(slice1, NumLevels, myPerm.CurrentParticles,0);
//       }
//       else 
// 	acceptBisect=false;
//       if (acceptBisect){
// 	//	cerr<<"I've accepted!"<<endl;
// 	PathData.AcceptMove(slice1,slice2,myPerm.CurrentParticles);
//       }
//       else 
// 	PathData.RejectMove(slice1, slice2,myPerm.CurrentParticles);
//     }
//   }    
//   void Read(IOSectionClass &in)
//   {
//     assert(in.ReadVar("name",Name));
//     assert(in.ReadVar("NumLevels",NumLevels));
//   }
   PermMove(PathDataClass &myPathData, IOSectionClass outSection) : 
     ParticleMoveClass (myPathData, outSection)
   {
     /* do nothing for now */
   }


};





#endif
