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

#include "DynamicStructureFactor.h"
#include <utility>
// #include <multimap>

void DynamicStructureFactorClass::Read(IOSectionClass& in)
{
  cerr<<"READING DYNAMIC STRUCTURE FACTOR"<<endl;
  
  ObservableClass::Read(in);
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

  TotalCounts=0;

  

  //These are kVecs that wouldn't be calculated given the kcutoff
  Array<double,2> tempkvecs;
  Array<double,1> kMagRange;
  if (in.ReadVar("AdditionalkVecs",tempkvecs)){
    assert(tempkvecs.extent(1)==NDIM);
    Additionalkvecs.resize(tempkvecs.extent(0));
    for (int kvec=0;kvec<tempkvecs.extent(0);kvec++)
      for (int dim=0;dim<NDIM;dim++){
	Additionalkvecs(kvec)[dim]=tempkvecs(kvec,dim);
      }
  }
  else if (in.ReadVar("kMagRange",kMagRange)){
    vector<dVec> tempkvecs_vec;
    assert(kMagRange.size()==2);
#if NDIM==2
    int maxI=int(kMagRange(1)*PathData.Path.GetBox()[0]/(2.0*M_PI))+1;
    int maxJ=int(kMagRange(1)*PathData.Path.GetBox()[1]/(2.0*M_PI))+1;
    cerr<<"Maxes are "<<maxI<<" "<<maxJ<<endl;
    for (int i=0;i<maxI;i++)
      for (int j=0;j<maxJ;j++){
	double kmag=sqrt(
			 (2*i*M_PI/PathData.Path.GetBox()[0])*
			 (2*i*M_PI/PathData.Path.GetBox()[0])+
			 (2*j*M_PI/PathData.Path.GetBox()[1])*
			 (2*j*M_PI/PathData.Path.GetBox()[1]));
	//	cerr<<kMagRange(0)<<" "<<kmag<<" "<<kMagRange(1)<<" "<<i<<" "<<j<<endl;
	if (kMagRange(0)<kmag && kmag<kMagRange(1)){
	  dVec kVec((2*i*M_PI/PathData.Path.GetBox()[0]),(2*j*M_PI/PathData.Path.GetBox()[1]));
	  tempkvecs_vec.push_back(kVec);
	  kVec[0]=-kVec[0];
	  tempkvecs_vec.push_back(kVec);
	}
      }
    Additionalkvecs.resize(tempkvecs_vec.size());
    //    cerr<<"Additional kvecs are "<<endl;
    for (int kvec=0;kvec<tempkvecs_vec.size();kvec++){
      Additionalkvecs(kvec)=tempkvecs_vec[kvec];
      //      cerr<<Additionalkvecs(kvec)[0]<<" "<<Additionalkvecs(kvec)[1]<<endl;
    }
    cerr<<"done"<<endl;
#else
	cerr<<"DOES NOT SUPPORT 3d yet"<<endl;
#endif
  }	  
  else{
    Additionalkvecs.resize(0);
  }
  AdditionalRho_k.resize(PathData.Path.NumTimeSlices(),1, Additionalkvecs.size()); 
  

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
  Sk.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
    Fkt.resize(PathData.Path.TotalNumSlices,PathData.Path.kVecs.size());
  rho_k_real.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
  rho_k_imag.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
  Sk=0;
  cerr<<"Done READING DYNAMIC STRUCTURE FACTOR"<<endl;
}



void DynamicStructureFactorClass::WriteInfo()
{
  PathClass &Path= PathData.Path;
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  ObservableClass::WriteInfo();
  Array<double,2> kVecArray(kVecs.size()+Additionalkvecs.size(),NDIM);
  Array<double,1> kMagArray(kVecs.size()+Additionalkvecs.size());

  for (int ki=0; ki < kVecs.size(); ki++) {
    dVec &k = kVecs(ki);
    kMagArray(ki) = sqrt(dot(k,k));
    for (int j=0; j<NDIM; j++)
      kVecArray(ki,j) = k[j];
  }
  for (int ki=kVecs.size();ki<kVecs.size()+Additionalkvecs.size();ki++){
    dVec &k = Additionalkvecs(ki-kVecs.size());
    kMagArray(ki) = sqrt(dot(k,k));
    for (int j=0; j<NDIM; j++)
      kVecArray(ki,j) = k[j];
  }
  IOSection.WriteVar("kVecs", kVecArray);
  ///We now accumulate the structure factor one at a time
  IOSection.WriteVar("Cumulative","False");
  /// Output data for plotting in analysis code
  IOSection.WriteVar("x", kMagArray);
  IOSection.WriteVar("xlabel", "|k|");
  IOSection.WriteVar("ylabel", "S(k)");
  IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
  IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
  IOSection.WriteVar("Type","CorrelationFunction");
}


void DynamicStructureFactorClass::WriteBlock()
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  Array<double,2> SkSum(PathData.Path.TotalNumSlices,kVecs.size());
  SkSum=0.0;
  double norm=0.0;
  int num1 = PathData.Path.Species(Species1).NumParticles;
  int num2 = PathData.Path.Species(Species1).NumParticles;
  norm = PathData.Path.TotalNumSlices*TotalCounts * sqrt((double)num1*num2);
  PathData.Path.Communicator.Sum(Fkt, SkSum);
  if (PathData.Path.Communicator.MyProc()==0) 
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
  Array<double,2> SofkArray(SkSum.extent(0),SkSum.extent(1));
  for (int ki=0; ki<kVecs.size(); ki++)
    for (int shift=0;shift<PathData.Path.TotalNumSlices;shift++)
      SofkArray(shift,ki) = (double) Fkt(shift,ki) / norm;
  //  for (int ki=kVecs.size();ki<kVecs.size()+Additionalkvecs.size();ki++)
  //    SofkArray(Range::all(),ki) = (double) SkSum(ki) / norm;
  SofkVar.Write(SofkArray);
  ///Clear the structure factor counts
  Fkt=0;
  Sk=0;
  rho_k_real=0;
  rho_k_imag=0;
  TotalCounts=0;
  SkMax=0;
  MaxkVec=0;
}




void DynamicStructureFactorClass::Accumulate()
{
  cerr<<"ACCUMULATE DYNAMIC STRUCTURE FACTOR"<<endl;
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  //  cerr<<"I have been told to accumulate"<<endl;
  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);


  if (!PathData.Path.LongRange) {
    for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
      PathData.Path.CalcRho_ks_Fast(slice, Species1);
    if (Species2 != Species1)
      for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
	PathData.Path.CalcRho_ks_Fast(slice, Species2);
  }
  if (Additionalkvecs.extent(0)!=0){
    assert(Species1==Species2);
    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++)
      PathData.Path.CalcRho_ks_Slow(slice,Species1,
				    Additionalkvecs,
				    AdditionalRho_k);
  }

  //  int totalSlices=Path.TotalNumSlices;
  //  Array<double,3> RhokArray;
  //  RhokArray.resize(Path.kIndices.size(),totalSlices,NDIM);
  cerr<<"INTO LOOP"<<endl;
  TotalCounts++;
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    for (int shift=0;shift<PathData.Path.TotalNumSlices;shift++){
      for (int ki=0; ki<kVecs.size(); ki++) {
	double a = PathData.Path.Rho_k(slice, Species1, ki).real();
	double b = PathData.Path.Rho_k(slice, Species1, ki).imag();
	double c = PathData.Path.Rho_k((slice+shift) % (PathData.Path.NumTimeSlices()-1), Species2, ki).real();
	double d = PathData.Path.Rho_k((slice+shift) % (PathData.Path.NumTimeSlices()-1), Species2, ki).imag();
	double sk=a*c+b*d;
	Fkt(shift,ki)+=sk;
      }
    }
  }
  cerr<<"Done ACCUMULATE DYNAMIC STRUCTURE FACTOR"<<endl;
}
void DynamicStructureFactorClass::Clear()
{
  Fkt=0;
  TotalCounts=0;
}


void DynamicStructureFactorClass::Initialize()
{
  Fkt=0;
  TotalCounts = 0;
  TimesCalled=0;
  SkMax=0;
  MaxkVec=0;
  

}


