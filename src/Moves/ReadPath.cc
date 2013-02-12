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

#include "ReadPath.h"
///BUG: Doesn't load permutations yet
void 
ReadPathClass::MakeMove()
{
  Array<double,3> PathVec;
  IOVar->Read(PathVec,currConfig,Range::all(),Range::all(),Range::all());
  cerr<<PathVec.extent(0)<<" "<<PathVec.extent(1)<<" "<<PathVec.extent(2)<<endl;
  cerr<<PathData.Path.NumTimeSlices()<<endl;
  assert(PathVec.extent(0)==PathData.Path.NumParticles());
  assert(PathVec.extent(1)==PathData.Path.NumTimeSlices()-1);
  assert(PathVec.extent(2)==NDIM);
  
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
      dVec pos;
      for (int dim=0;dim<NDIM;dim++)
 	pos(dim)=PathVec(ptcl,slice,dim);
      PathData.Path.SetPos(slice,ptcl,pos);
      //   cerr<<PathVec.extent(0)<<" "<<PathVec.extent(1)<<" "<<PathVec.extent(2)<<endl;
    }
  }
  cerr<<"Going to accept"<<endl;
  PathData.AcceptMove(0,PathData.Path.NumTimeSlices()-1,ActiveParticles);
  currConfig++;
  if (currConfig>=NumConfigs){
    cerr<<"You are trying to read configurations that aren't there"<<endl;
    assert(1==2);
  }
}

void
ReadPathClass::Read(IOSectionClass &input)
{
  cerr<<"I'm being read"<<endl;
  string fileName;
  assert(input.ReadVar("File",fileName));
  IOSectionClass in;
  assert (in.OpenFile(fileName.c_str()));
  assert(in.OpenSection("Observables"));
  assert(in.OpenSection("PathDump"));

  cerr<<"This far"<<endl;
  IOVar = in.GetVarPtr("Path");
  Array<double,1> checkSize;
  IOVar->Read(checkSize,Range::all(),0,0,0);
  NumConfigs=checkSize.size();
  currConfig=0;
  ActiveParticles.resize(PathData.Path.NumParticles());
  for (int i=0;i<ActiveParticles.size();i++)
    ActiveParticles(i)=i;
}


