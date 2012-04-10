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

#include "DistanceToHead.h"


// Fix to include final link between link M and 0
void HeadLocClass::Accumulate()
{
  
  int closestHeadLoc=0;
  int closestTailLoc=0;
  ///Multiplication by 5 just to make sure closestHead and closestTail
  ///are larger then you could ever have the distance to the head or
  ///tail from any fixed location
  double closestHead=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
  double closestTail=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
  dVec headDisp,tailDisp;
  double headDist2,tailDist2; 
  int openSlice=PathData.Path.OpenLink;
  int headPtcl=PathData.Path.OpenPtcl;
  ///In non-open loops mode this particle wouldn't exist.  In open
  ///loop mood it stores the location of the tail.
  int tailPtcl=PathData.Path.NumParticles(); 				   
  dVec headLoc=PathData.Path(openSlice,headPtcl);
  dVec tailLoc=PathData.Path(openSlice,tailPtcl);
  for (int counter=0;counter<FixedLoc.size();counter++){
    headDisp=FixedLoc(counter)-headLoc;
    tailDisp=FixedLoc(counter)-tailLoc;
    PathData.Path.PutInBox(headDisp);
    PathData.Path.PutInBox(tailDisp);
    headDist2=dot(headDisp,headDisp);
    tailDist2=dot(tailDisp,tailDisp);
    if (headDist2<closestHead){
      closestHead=headDist2;
      closestHeadLoc=counter;
    }
    if (tailDist2<closestTail){
      closestTail=tailDist2;
      closestTailLoc=counter;
    }
  }
  HeadLoc(closestHeadLoc)++;
  TailLoc(closestTailLoc)++;
  NumSamples++;
}

void HeadLocClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  Array<int,1> headLocSum(HeadLoc.size());
  headLocSum=0;
  Array<int,1> tailLocSum(TailLoc.size());
  tailLocSum=0;
  PathData.Path.Communicator.Sum(HeadLoc,headLocSum);
  PathData.Path.Communicator.Sum(TailLoc,tailLocSum);
  Array<double,1> headLocDouble(HeadLoc.size());
  headLocDouble=0.0;
  Array<double,1> tailLocDouble(TailLoc.size());
  tailLocDouble=0.0;
  for (int counter=0;counter<headLocSum.size();counter++){
    headLocDouble(counter)=headLocSum(counter)*norm;
    tailLocDouble(counter)=tailLocSum(counter)*norm;
  }
  HeadLocVar.Write(headLocDouble);
  TailLocVar.Write(tailLocDouble);
  HeadLoc=0;
  TailLoc=0;
  NumSamples = 0;
}

void HeadLocClass::Read(IOSectionClass &in)
{  
  int numFixedPoints;
  ObservableClass::Read(in);
  assert(in.ReadVar("NumFixedPoints",numFixedPoints));
  HeadLoc.resize(numFixedPoints);
  TailLoc.resize(numFixedPoints);
  FixedLoc.resize(numFixedPoints);
  HeadLoc=0;
  TailLoc=0;
  Array<double,2> positions;
  assert(in.ReadVar("LocationsToCompare",positions));
  ///Verify you used the right number of points to compare against
  assert(positions.extent(0)==HeadLoc.size());
  assert(positions.extent(1)==NDIM);
  dVec pos;
  for (int loc=0;loc<FixedLoc.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    FixedLoc(loc) = pos;
  }      
  
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}



void HeadLocClass::WriteInfo()
{


}
