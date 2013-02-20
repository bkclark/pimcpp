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

#include "StructureFactor.h"
#include <utility>
#include <map>

using namespace std;

void StructureFactorClass::Read(IOSectionClass& in)
{
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
      for (int j=0;j<maxJ;j++) {
        double kmag=sqrt((2*i*M_PI/PathData.Path.GetBox()[0])*
                         (2*i*M_PI/PathData.Path.GetBox()[0])+
                         (2*j*M_PI/PathData.Path.GetBox()[1])*
                         (2*j*M_PI/PathData.Path.GetBox()[1]));
        if (kMagRange(0)<kmag && kmag<kMagRange(1)) {
          dVec kVec((2*i*M_PI/PathData.Path.GetBox()[0]),(2*j*M_PI/PathData.Path.GetBox()[1]));
          tempkvecs_vec.push_back(kVec);
          kVec[0]=-kVec[0];
          tempkvecs_vec.push_back(kVec);
        }
      }
    Additionalkvecs.resize(tempkvecs_vec.size());
    for (int kvec=0;kvec<tempkvecs_vec.size();kvec++) {
      Additionalkvecs(kvec)=tempkvecs_vec[kvec];
    }
    cerr<<"done"<<endl;
#else
    cerr<<"DOES NOT SUPPORT 3d yet"<<endl;
#endif
  } else
    Additionalkvecs.resize(0);
  AdditionalRho_k.resize(PathData.Path.NumTimeSlices(),1, Additionalkvecs.size());

  /// Setup k Vecs and RhoK
  PathData.Path.SetupkVecs(in);

  // Sign Tracking
  if(!in.ReadVar("TrackSign", TrackSign))
    TrackSign = 0;

  Sk.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
  rho_k_real.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
  rho_k_imag.resize(PathData.Path.kVecs.size()+Additionalkvecs.size());
  Sk=0;
}


void StructureFactorClass::WriteInfo()
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


void StructureFactorClass::WriteBlock()
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  Array<double,1> SkSum(kVecs.size()+Additionalkvecs.size());
  SkSum=0.0;
  double norm=0.0;
  int num1 = PathData.Path.Species(Species1).NumParticles;
  int num2 = PathData.Path.Species(Species1).NumParticles;
  norm = PathData.Path.TotalNumSlices*TotalCounts * sqrt((double)num1*num2);
  SkMaxVar.Write(SkMax);

  Array<double,1> rho_k_realSum(kVecs.size()+Additionalkvecs.size());
  rho_k_realSum=0.0;
  PathData.Path.Communicator.Sum(rho_k_real, rho_k_realSum);
  rho_k_realSum=rho_k_realSum/norm;
  rho_k_realVar.Write(rho_k_realSum);

  rho_k_realSum=0.0;
  PathData.Path.Communicator.Sum(rho_k_imag, rho_k_realSum);
  rho_k_realSum=rho_k_realSum/norm;
  rho_k_imagVar.Write(rho_k_realSum);

    

  PathData.Path.Communicator.Sum(Sk, SkSum);
  if (PathData.Path.Communicator.MyProc()==0) 
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
    }
  Array<double,1> SofkArray(SkSum.size());
  for (int ki=0; ki<kVecs.size(); ki++)
    SofkArray(ki) = (double) SkSum(ki) / norm;
  for (int ki=kVecs.size();ki<kVecs.size()+Additionalkvecs.size();ki++)
    SofkArray(ki) = (double) SkSum(ki) / norm;
  SofkVar.Write(SofkArray);
  ///Clear the structure factor counts
  Sk=0;
  rho_k_real=0;
  rho_k_imag=0;
  TotalCounts=0;
  SkMax=0;
  MaxkVec=0;
}




void StructureFactorClass::Accumulate()
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  //  cerr<<"I have been told to accumulate"<<endl;
  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);

  TotalCounts++;

  double FullWeight = CalcFullWeight();

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
      PathData.Path.CalcRho_ks_Slow(slice,Species1,Additionalkvecs,AdditionalRho_k);
  }
  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    multimap<double,double > kList;
    for (int ki=0; ki<kVecs.size(); ki++) {
      double a = PathData.Path.Rho_k(slice, Species1, ki).real();
      double b = PathData.Path.Rho_k(slice, Species1, ki).imag();
      double c = PathData.Path.Rho_k(slice, Species2, ki).real();
      double d = PathData.Path.Rho_k(slice, Species2, ki).imag();
      // \f$ Sk(ki) :=  Sk(ki) + \Re(rho^1_k * rho^2_{-k}) \f
      double sk=a*c+b*d;
      rho_k_real(ki) += FullWeight*a;
      rho_k_imag(ki) += FullWeight*b;
      if (sk>SkMax){
       SkMax=sk;
       MaxkVec=kVecs(ki);
      }
      double kMag=sqrt(kVecs(ki)[0]*kVecs(ki)[0]+kVecs(ki)[1]*kVecs(ki)[1]);
      kList.insert(pair<double,double> (kMag,sk));
                  //      cerr<<slice<<" "<<ki<<" "<<sk<<endl;
      Sk(ki) += FullWeight*sk;
    }
    for (int ki=kVecs.size();ki<kVecs.size()+Additionalkvecs.size();ki++){
      int kk=ki-kVecs.size();
      double a = (AdditionalRho_k(slice, Species1, kk)).real();
      double b = AdditionalRho_k(slice, Species1, kk).imag();
      double c = AdditionalRho_k(slice, Species2, kk).real();
      double d = AdditionalRho_k(slice, Species2, kk).imag();

      // \f$ Sk(ki) :=  Sk(ki) + \Re(rho^1_k * rho^2_{-k}) \f
      double sk=a*c+b*d;
      Sk(ki) += FullWeight*sk;
      rho_k_real(ki) += FullWeight*a;
      rho_k_imag(ki) += FullWeight*b;

    }
  }

}


void StructureFactorClass::Clear()
{
  Sk=0;
  rho_k_real=0;
  rho_k_imag=0;

  TotalCounts=0;
}


void StructureFactorClass::Calculate()
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

}
void StructureFactorClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;
  SkMax=0;
  MaxkVec=0;

}


