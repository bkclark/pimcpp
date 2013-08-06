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

#include "CenterofMassMove.h"


void CenterOfMassMoveClass::MakeMove()
{
  PathClass &Path = PathData.Path;
  //  cerr<<"The current center of mass is "<<Path.CenterOfMass;
  dVec toShift=-1* (Path.CenterOfMass/((Path.NumTimeSlices()-1)*Path.NumParticles()));
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
      Path(slice,ptcl)=Path(slice,ptcl)+toShift;
    }
  int slice1=0;
  int slice2=PathData.Path.NumTimeSlices()-1;
  
  PathData.AcceptMove(slice1,slice2,ActiveParticles);
  //  cerr<<"The post center of mass is "<<Path.CenterOfMass;
  return;

  PathData.AcceptMove(slice1,slice2,ActiveParticles);


  dVec center_of_mass=0.0;
  dVec zero=0.0;
  int totalCount=0;
  //Calculate center of mass
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    dVec sliceZero=Path(0,ptcl);
    dVec pathCenter=0.0;
    for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
      pathCenter+=PathData.Path.Velocity(0,slice,ptcl);
    }
    pathCenter=pathCenter/(PathData.Path.NumTimeSlices()-2);
    dVec tempDist=sliceZero+pathCenter;
    center_of_mass+=PathData.Path.MinImageDisp(zero,tempDist);
    totalCount++;
  }
  //  center_of_mass=PathData.Path.MinImageDisp(center_of_mass,zero);
  //  cerr<<"My pre-center of mass is "<<center_of_mass<<endl;
  center_of_mass=center_of_mass/(totalCount+0.0);
  if (firstTime){
    original_center_of_mass=center_of_mass;
    firstTime=false;
  }
  center_of_mass=center_of_mass-original_center_of_mass;
  if (abs(center_of_mass[0])>1.0 || abs(center_of_mass[1])>1.0 ||
      abs(center_of_mass[2])>1.0){
      
  cerr<<"My center of mass is "<<center_of_mass<<endl;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
      Path(slice,ptcl)=Path(slice,ptcl)-(1.0/5.0)*center_of_mass;
      //      dVec cm=center_of_mass/5.0;
      //      cm=-cm;
      //      Path(slice,ptcl)=Path(slice,ptcl)+cm;
    }
  
  int slice1=0;
  int slice2=PathData.Path.NumTimeSlices()-1;
  
  PathData.AcceptMove(slice1,slice2,ActiveParticles);

	  

  center_of_mass=0.0;
  zero=0.0;
  totalCount=0;

  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    dVec sliceZero=Path(0,ptcl);
    dVec pathCenter=0.0;
    for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
      pathCenter+=PathData.Path.Velocity(0,slice,ptcl);
    }
    pathCenter=pathCenter/(PathData.Path.NumTimeSlices()-2);
    dVec tempDist=sliceZero+pathCenter;
    center_of_mass+=PathData.Path.MinImageDisp(zero,tempDist);
    totalCount++;
  }
  //  center_of_mass=PathData.Path.MinImageDisp(center_of_mass,zero);
  //  cerr<<"My pre-center of mass is "<<center_of_mass<<endl;
  center_of_mass=center_of_mass/(totalCount+0.0);
  if (firstTime){
    original_center_of_mass=center_of_mass;
    firstTime=false;
  }
  center_of_mass=center_of_mass-original_center_of_mass;



  cerr<<"My post center of mass is "<<center_of_mass<<endl;
  }

}

void CenterOfMassMoveClass::Read(IOSectionClass &in)
{
  string typeCheck;
  assert(in.ReadVar("Type",typeCheck));
  assert(typeCheck=="CenterOfMass");
  assert(in.ReadVar("Name",Name));
  ActiveParticles.resize(PathData.Path.NumParticles());
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
    ActiveParticles(ptcl)=ptcl;
  firstTime=true;
}
