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


#include "DavidLongRangeClassYk.h"
#include "../PathDataClass.h"
#include <set>
#include "sys/time.h"


bool DavidLongRangeClassYk::fequals(double a,double b,double tol)
{
  return abs(a-b)<tol;
}

bool DavidLongRangeClassYk::vecEquals(dVec &a, dVec &b,double tol)
{
  bool equals=true;
  for (int dim=0;dim<NDIM;dim++)
    equals= equals && fequals(a[dim],b[dim],tol);
  return equals;
}

void DavidLongRangeClassYk::Build_MultipleSpecies()
{
  // Volume
  double vol=1.0;
  for (int dim=0; dim<NDIM; dim++)
    vol*=Path.GetBox()[dim];

  for (int i=0; i<duk.extent(0); i++)
    for (int j=0; j<duk.extent(1); j++)
      duk(i,j)=100.0;

  for (int iPair=0; iPair<PairArray.size(); iPair++) {
    DavidPAClass &pa(*(DavidPAClass*)PairArray(iPair));
    for (int i=0; i<pa.kVals.size(); i++){
      double k=pa.kVals(i);
      if (fequals(0.0,k,1e-10))
        yk_zero(iPair) = pa.uk_long(i)/vol;
      bool found=false;
      for (int j=0; j<Path.kVecs.size(); j++) {
        if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-4)) {
          Vlong_k(iPair,j) = pa.uk_long(i)/vol;
          found=true;

          if (pa.HasLongRange) {
            uk(iPair,j) = Vlong_k(iPair,j)*Path.tau;
            duk(iPair,j) = Vlong_k(iPair,j);
          } else {
            uk(iPair,j) = 0.;
            duk(iPair,j) = 0.;
          }
        }
      }
    }
  }
}

void DavidLongRangeClassYk::BuildRPA_MultipleSpecies()
{
  if (PathData.Path.Communicator.MyProc() == 0)
    cout << Path.CloneStr << " Using Long Range RPA correction " << endl;
  for (int speciesNum = 0; speciesNum < PairArray.size(); speciesNum++) {
    for (int i=0;i<uk.extent(1);i++)
      uk(speciesNum,i)=-1;
    double vol=1.0;
    for (int dim=0;dim<NDIM;dim++)
      vol*=Path.GetBox()[dim];
    DavidPAClass &pa(*(DavidPAClass*)PairArray(speciesNum));
    int ncomps=Path.Species(speciesNum).NumParticles;
    double lambda=Path.Species(speciesNum).lambda;
    for (int i=0;i<pa.kVals.size();i++){
      double k=pa.kVals(i);
      if (fequals(0.0,k,1e-10))
        yk_zero(speciesNum)=pa.uk_long(i)/(vol);
    }
    for (int i=0;i<pa.kVals.size();i++){
      double k=pa.kVals(i);
      double arg=1.0+ncomps*2.0*pa.uk_long(i)/(lambda*k*k*vol);
      double theta=0.5*k*k*lambda*Path.tau;
      double q;
      double s;
      double tn;
      if (arg<0){
        q=sqrt(-arg);
        tn=tan(q*theta);
        s=q*(1-q*tn)/(q+tn);
      }
      else if (arg==0)
        s=1.0/(1.0+theta);
      else {
        q=sqrt(arg);
        tn=tanh(theta*q);
        s=q*(1.0+q*tn)/(q+tn);
      }
      for (int j=0;j<Path.kVecs.size();j++){
        if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-10)){
           uk(speciesNum, j)=(-1.0+s)/ncomps;
           //cout << Path.kVecs(j) << " " << uk(speciesNum,j) << " " << Path.tau*pa.uk_long(i)/(vol) << endl;
           duk(speciesNum, j)=pa.uk_long(i)/(vol)-lambda*k*k*uk(speciesNum, j)*((1.0+0.5*ncomps*uk(speciesNum, j)));
           Vlong_k(speciesNum, j)=pa.uk_long(i)/(vol);
           //cerr<<"VLONG IS "<<Vlong_k(speciesNum, j)<<" "<<pa.uk_long(i)<<" "<<speciesNum<<endl;
        }
      }
    }
    cout<<PathData.Path.Communicator.MyProc()<<" RPA Built: species("<<speciesNum<<") lambda("<<lambda<<") N("<<ncomps<<") Vol("<<vol<<")"<<endl;
    //for (int i=0;i<uk.extent(1);i++)
    //  cerr<<"KVecs: "<<sqrt(blitz::dot(Path.kVecs(i),Path.kVecs(i)))<<" "<<uk(speciesNum, i)<<" "<<duk(speciesNum, i)<<endl;
  }
}

void DavidLongRangeClassYk::BuildRPA_SingleSpecies()
{
  for (int i=0;i<uk.extent(1);i++)
    uk(0,i)=-1;
  int speciesNum=0;
  double vol=1.0;
  for (int dim=0;dim<NDIM;dim++)
    vol*=Path.GetBox()[dim];
  DavidPAClass &pa(*(DavidPAClass*)PairArray(speciesNum));
  int ncomps=Path.Species(speciesNum).NumParticles;
  double lambda=Path.Species(speciesNum).lambda;
  //cerr<<"lambda is "<<lambda<<endl;
  //cerr<<"ncomps is "<<ncomps<<endl;
  //cerr<<"Vol is "<<vol<<endl;
  for (int i=0;i<pa.kVals.size();i++){
    double k=pa.kVals(i);
    if (fequals(0.0,k,1e-10))
      yk_zero(0)=pa.uk_long(i)/(vol);
  }
  for (int i=0;i<pa.kVals.size();i++){
    double k=pa.kVals(i);
    double arg=1.0+ncomps*2.0*pa.uk_long(i)/(lambda*k*k*vol);
    double theta=0.5*k*k*lambda*Path.tau;
    double q;
    double s;
    double tn;
    if (arg<0){
      q=sqrt(-arg);
      tn=tan(q*theta);
      s=q*(1-q*tn)/(q+tn);
    }
    else if (arg==0)
      s=1.0/(1.0+theta);
    else {
      q=sqrt(arg);
      tn=tanh(theta*q);
      s=q*(1.0+q*tn)/(q+tn);
    }
    for (int j=0;j<Path.kVecs.size();j++){
      if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-10)){
        uk(speciesNum, j)=(-1.0+s)/ncomps;
        duk(speciesNum, j)=pa.uk_long(i)/(vol)-lambda*k*k*uk(speciesNum, j)*((1.0+0.5*ncomps*uk(speciesNum, j)));
        Vlong_k(speciesNum, j)=pa.uk_long(i)/(vol);
        //cerr<<"VLONG IS "<<Vlong_k(speciesNum, j)<<endl;
      }
    }
  }
  //cerr<<"I have built the rpa"<<endl;
  //for (int i=0;i<uk.extent(1);i++)
  //  cerr<<"KVecs: "<<sqrt(blitz::dot(Path.kVecs(i),Path.kVecs(i)))<<" "<<uk(speciesNum, i)<<" "<<duk(speciesNum, i)<<endl;
}

void DavidLongRangeClassYk::WriteInfo(IOSectionClass &out)
{
  out.WriteVar ("Type", "David_Long_Range");
}

void DavidLongRangeClassYk::Read(IOSectionClass &in)
{
  // Set flags
  Path.LongRange = true;
  Path.DavidLongRange = true;

  // Setup k Vecs and RhoK
  Path.SetupkVecs(in);

  TimeSpent=0;
  for (int pai=0;pai<PairArray.size();pai++){
    DavidPAClass &pa(*((DavidPAClass*)PairArray(pai)));
    //cerr<<"the long range is "<<pa.LongRangeDim<<endl;
    //    assert(pa.LongRangeDim==NDIM);
    //for (int dim=0;dim<NDIM;dim++){
    //  cerr<<pa.LongRangeBox(dim)<<" "<<Path.GetBox()[dim]<<endl;
    //  //      assert(fabs(pa.LongRangeBox(dim)-Path.GetBox()[dim])<1e-5);
    //}
    int specNum1=0;
    while (Path.Species(specNum1).Type!=pa.Particle1.Name){
      specNum1++;
      assert(specNum1<Path.NumSpecies());
    }
    int specNum2=0;
    while (Path.Species(specNum2).Type!=pa.Particle2.Name){
      specNum2++;
      assert(specNum2<Path.NumSpecies());
    }
    //cerr<<"MASSES: "<<pa.LongRangeMass1<<" "<<pa.LongRangeMass2<<" "<<Path.Species(specNum1).lambda<<" "<<Path.Species(specNum2).lambda<<endl;
    //    assert(fabs(pa.LongRangeMass1-Path.Species(specNum1).lambda)<1e-10);
    //    assert(fabs(pa.LongRangeMass2-Path.Species(specNum2).lambda)<1e-10);
  }

  uk.resize(PairArray.size(), Path.kVecs.size());
  duk.resize(PairArray.size(), Path.kVecs.size());
  Vlong_k.resize(PairArray.size(), Path.kVecs.size());
  yk_zero.resize(PairArray.size());
  Spec2Index.resize(PairMatrix.extent(0),PairMatrix.extent(1));

  // Check to see if doing RPA
  if (!in.ReadVar ("UseRPA", UseRPA))
    UseRPA = false;
  if (PathData.Actions.UseRPA) {
    BuildRPA_MultipleSpecies();
  } else {
    Build_MultipleSpecies();
  }
}

/// Calculates the long range part of the action using David's breakup.
/// The short range part must be supplied as a dm file without the long
/// range part in it.  It ignores active particles.
double DavidLongRangeClassYk::SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  int skip = 1<<level;

  set<int> speciesList;
  for(int p=0; p<activeParticles.size(); p++) {
    int ptcl = activeParticles(p);
    int spec = Path.ParticleSpeciesNum(ptcl);
    speciesList.insert(speciesList.begin(), spec);
  }

  bool only_do_inclusive=false;
  if (slice1==0 && slice2==PathData.Path.NumTimeSlices()-1)
    only_do_inclusive = false;
  else
    only_do_inclusive = true;
  if (GetMode() == NEWMODE) {
    if (!only_do_inclusive)
      Path.UpdateRho_ks(slice1, slice2-skip, activeParticles, level);
    else
      Path.UpdateRho_ks(slice1+skip, slice2-skip, activeParticles, level);
  }

  int startSlice = slice1;
  int endSlice = slice2-skip;
  if (only_do_inclusive) {
    startSlice = slice1+skip;
    endSlice = slice2-skip;
  }

  double total = 0;
  double factor = 1.0;
  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=startSlice; slice<=endSlice; slice+=skip) {
      for(set<int>::iterator it = speciesList.begin(); it!=speciesList.end(); it++) {
        int species = *it;
        double rhok2 = mag2(Path.Rho_k(slice,species,ki));
        total +=  factor*rhok2*uk(PairIndex(species,species),ki);
      }
    }
  }

  // Cross-terms for Multiple Species
  factor = 2.0;
  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=startSlice; slice<=endSlice; slice+=skip) {
      for (int species0=0; species0<Path.NumSpecies()-1; species0++) {
        for (int species1=species0+1; species1<Path.NumSpecies(); species1++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          total += factor*rhok2*uk(PairIndex(species0,species1),ki);
        }
      }
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return total;
}

double DavidLongRangeClassYk::d_dBeta (int slice1, int slice2,  int level)
{
  double total=0.0;
  double factor=1.0;

  for (int slice=slice1; slice<=slice2; slice++) {
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
         double rhok2 = mag2(Path.Rho_k(slice,species,ki));
         sliceTotal +=  factor*rhok2*duk(PairIndex(species,species),ki);
      }
    }
    total += sliceTotal;
  }

  // Cross-terms for Multiple Species
  for (int slice=slice1; slice<=slice2; slice++) {
    double sliceTotal=0.0;
    if ((slice==slice1) || (slice==slice2))
      factor = 1.0;
    else
      factor = 2.0;
    for (int species0=0; species0<Path.NumSpecies()-1; species0++) {
      Path.CalcRho_ks_Fast(slice,species0);
      for (int species1=species0+1; species1<Path.NumSpecies(); species1++) {
        Path.CalcRho_ks_Fast(slice,species1);
        for (int ki = 0; ki < Path.kVecs.size(); ki++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          sliceTotal += factor*rhok2*duk(PairIndex(species0,species1),ki);
        }
      }
    }
    total += sliceTotal;
  }

  return total;

}

///Not really d_dbeta but total energy
double DavidLongRangeClassYk::V (int slice1, int slice2,  int level)
{
  double total=0.0;
  double factor=1.0;
  for (int slice=slice1; slice<=slice2; slice++) {
    double sliceTotal = 0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
         double rhok2 = mag2(Path.Rho_k(slice,species,ki));
         sliceTotal += factor*rhok2*Vlong_k(PairIndex(species,species),ki);
      }
    }
    total += sliceTotal;
  }

  // Cross-terms for Multiple Species
  for (int slice=slice1;slice<=slice2;slice++) {
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 1.0;
    else
      factor = 2.0;
    for (int species0 = 0; species0 < Path.NumSpecies()-1; species0++) {
      Path.CalcRho_ks_Fast(slice,species0);
      for (int species1 = species0+1; species1 < Path.NumSpecies(); species1++) {
        Path.CalcRho_ks_Fast(slice,species1);
        for (int ki = 0; ki < Path.kVecs.size(); ki++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          total += factor*rhok2*Vlong_k(PairIndex(species0,species1),ki);
        }
      }
    }
  }

  return total;

}

string DavidLongRangeClassYk::GetName()
{
  return "DavidLongRange";
}
