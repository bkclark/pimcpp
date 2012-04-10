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

#include "OpenOrientation.h"


// Fix to include final link between link M and 0
void OpenOrientationClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }

  int openLink=PathData.Path.OpenLink;
  int openPtcl=PathData.Path.OpenPtcl;
  int numLinks=PathData.Path.NumTimeSlices()-1;
  dVec disp=0.0;
  int currSlice=openLink;
  int currPtcl=openPtcl;
  int nextSlice=-1;
  int nextPtcl=-1;
  

  //  cerr<<"hello"<<endl;
  while (nextSlice!=openLink || nextPtcl!=openPtcl){
    //    cerr<<nextSlice<<" "<<currSlice<<" "<<openLink<<endl;
    //    cerr<<nextPtcl<<" "<<currPtcl<<" "<<openPtcl<<endl;
    
    nextSlice = (currSlice + 1) % PathData.Path.NumTimeSlices();
    //    if (nextSlice==0)
    //      nextSlice=numLinks+1;
    if (currSlice==PathData.Join)
      nextPtcl=PathData.Path.Permutation(currPtcl);
    else 
      nextPtcl=currPtcl;
    dVec linkDisp=PathData.Path.VelocityBetweenPtcl(currSlice,currPtcl,nextSlice,nextPtcl);
    disp=disp+linkDisp;
    currSlice=nextSlice;
    currPtcl=nextPtcl;
  }
  double dist=sqrt(dot(disp,disp));
  R2+=disp(0)*disp(0)+disp(1)*disp(1);
  Z2+=disp(2);
  R2OverZ2+=R2/Z2;
  NumSamples++;
}

void OpenOrientationClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  R2Var.Write(R2*norm);
  Z2Var.Write(Z2*norm);
  R2OverZ2Var.Write(R2OverZ2*norm);
  R2=0.0;
  Z2=0.0;
  R2OverZ2=0.0;
  NumSamples = 0;
}

void OpenOrientationClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  R2=0.0;
  Z2=0.0;
  R2OverZ2=0.0;
  NumSamples=0;
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}



void OpenOrientationClass::WriteInfo()
{


}
