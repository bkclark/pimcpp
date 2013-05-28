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

#include "ParticleAverageLoc.h"

///This observable calculates the center of mass of all the particles
///in the. It currently will produce erroneous data if the particles
///in the system are permuted onto other particles (this is because it
///ignores the join and because it is unclear how to define the center
///of mass of such particles correctly.

void ParticleAverageLocClass::Accumulate()
{
  int firstPtcl=PathData.Path.Species(Species).FirstPtcl;
  int lastPtcl=PathData.Path.Species(Species).LastPtcl;

  dVec totalDispVec(0.0);
  dVec chainCtr;
  for (int ptcl=firstPtcl;ptcl<=lastPtcl;ptcl++){
    //We first get the vector to the first time slice of the particle
    dVec dispToPtcl=PathData.Path(0,ptcl);
    PathData.Path.PutInBox(dispToPtcl);
    chainCtr(0)=0.0;
    chainCtr(1)=0.0;
    chainCtr(2)=0.0;
    ///We now accumulate the center of mass with respect to this first particle
    for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
      dVec disp;
      disp=PathData.Path.Velocity(0,slice,ptcl);
      chainCtr+=disp;
    }
    chainCtr = chainCtr * (1.0/PathData.Path.NumTimeSlices());
    ///Now we add it onto the original vector
    chainCtr=chainCtr+dispToPtcl;
    PathData.Path.PutInBox(chainCtr);
    for (int dim=0;dim<NDIM;dim++)
      ParticleCenterOfMass(ptcl-firstPtcl,dim)+=chainCtr(dim);
  }
  NumSamples++;
}

  

void ParticleAverageLocClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  ParticleCenterOfMass *= norm;
  ParticleAverageLocVar.Write(ParticleCenterOfMass);
  ParticleCenterOfMass=0.0;
  NumSamples=0.0;
  ParticleAverageLocVar.Flush();
}

void ParticleAverageLocClass::Read(IOSectionClass &in)
{  
  string species1Name;
  Species=-1;
  assert(in.ReadVar("Species",species1Name));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species=spec;
    }
  }
  assert(Species!=-1);
  int numParticles=
    PathData.Path.Species(Species).LastPtcl-PathData.Path.Species(Species).FirstPtcl+1;
  
  ParticleCenterOfMass.resize(numParticles,NDIM);

  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }

  
  
 
}



void ParticleAverageLocClass::WriteInfo()
{
  

}
