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
#include "StructureReject.h"

///This has to be called after pathdata knows how many
///particles it has
void StructureRejectClass::Read(IOSectionClass& in)
{
  string species1Name;
  string species2Name;
  Species1=-1;
  Species2=-1;
  if (in.ReadVar("Species", species1Name)) {
    species2Name=species1Name;
  }
  else {
    assert(in.ReadVar("Species1",species1Name));
    assert(in.ReadVar("Species2",species2Name));
  }
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==species2Name){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);


  ///if it's not long range you haven't set the kvecs up yet and need to
  if (!PathData.Path.LongRange){//This is hackish..we use kcutoff
    ///to tell if you are long range and now we have to read 
    ///it in to get the structure factor corret.
    assert(in.ReadVar("kCutoff",PathData.Path.kCutoff));
  
#if NDIM==3    
    PathData.Path.SetupkVecs3D();
#endif
#if NDIM==2
    PathData.Path.SetupkVecs2D();
#endif
    PathData.Path.Rho_k.resize(PathData.Path.NumTimeSlices(), PathData.Path.NumSpecies(), PathData.Path.kVecs.size());
  }
  

  Sk.resize(PathData.Path.kVecs.size());
  Sk=0;


}

StructureRejectClass::StructureRejectClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double 
StructureRejectClass::SingleAction (int slice1, int slice2,
				    const Array<int,1> &changedParticles, 
				    int level)
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  TotalCounts++;
  if (!PathData.Path.LongRange) {
    for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
      PathData.Path.CalcRho_ks_Fast(slice, Species1);
    if (Species2 != Species1)
      for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
	PathData.Path.CalcRho_ks_Fast(slice, Species2);
  }

  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    for (int ki=0; ki<kVecs.size(); ki++) {
      double a = PathData.Path.Rho_k(slice, Species1, ki).real();
      double b = PathData.Path.Rho_k(slice, Species1, ki).imag();
      double c = PathData.Path.Rho_k(slice, Species2, ki).real();
      double d = PathData.Path.Rho_k(slice, Species2, ki).imag();
      // \f$ Sk(ki) :=  Sk(ki) + \Re(rho^1_k * rho^2_{-k}) \f
      Sk(ki) += a*c + b*d;	
    }
  }
  double max=Sk(0);
  for (int counter=0;counter<Sk.size();counter++){
    if (Sk(counter)<max)
      max=Sk(counter);
  }
  if (max>10.0){
    return 10e5;
  }
  else 
    return 0.0;

}



double StructureRejectClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}


string 
StructureRejectClass::GetName()
{
  return "StructureReject";
}
