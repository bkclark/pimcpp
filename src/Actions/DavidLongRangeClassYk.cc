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

#include "../MPI/Communication.h"
#include "DavidLongRangeClassYk.h"
#include "../PathDataClass.h"
#include <set>
#include  "sys/time.h"

bool fequals(double a,double b,double tol)
{
  return abs(a-b)<tol;
}

bool vecEquals(dVec &a, dVec &b,double tol)
{
  bool equals=true;
  for (int dim=0; dim<NDIM; dim++)
    equals = equals && fequals(a[dim],b[dim],tol);
  return equals;
}

void DavidLongRangeClassYk::Build_MultipleSpecies()
{
  //int speciesNum=0;
  double vol=1.0;
  for (int dim=0;dim<NDIM;dim++)
    vol*=Path.GetBox()[dim];
  for(int speciesNum=0; speciesNum<PairArray.size(); speciesNum++) {
    DavidPAClass &pa(*(DavidPAClass*)PairArray(speciesNum));
    //cerr<<"Vol is "<<vol<<endl;
    for (int i=0;i<pa.kVals.size();i++){
      double k=pa.kVals(i);
//       double arg=1.0+ncomps*2.0*pa.uk_long(i)/(lambda*k*k*vol);
//       double theta=0.5*k*k*lambda*Path.tau;
//       double q;
//       double s;
//       double tn;
//       if (arg<0){
//         q=sqrt(-arg);
//         tn=tan(q*theta);
//         s=q*(1-q*tn)/(q+tn);
//       }
//       else if (arg==0)
//         s=1.0/(1.0+theta);
//       else {
//         q=sqrt(arg);
//         tn=tanh(theta*q);
//         s=q*(1.0+q*tn)/(q+tn);
//       }
      for (int j=0;j<Path.kVecs.size();j++){
        if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-10)){
          Vlong_k(speciesNum, j)=pa.uk_long(i)/(vol);
          uk(speciesNum, j)=Vlong_k(speciesNum, j)*Path.tau;//(-1.0+s)/ncomps;
          duk(speciesNum, j)=Vlong_k(speciesNum, j);//pa.uk_long(i)/(vol)-lambda*k*k*uk(j)*((1.0+0.5*ncomps*uk(j)));
          //cout<<"VLONG IS "<<Vlong_k(j)<<endl;
        }
      }
    }
  }
}

void DavidLongRangeClassYk::BuildRPA_SingleType()
{
  int typeNum = 0;
  for (int speciesNum = 0; speciesNum < Path.NumSpecies(); speciesNum++) {
    for (int i=0;i<uk.extent(1);i++)
      uk(speciesNum,i)=-1;
    double vol=1.0;
    for (int dim=0;dim<NDIM;dim++)
      vol*=Path.GetBox()[dim];
    DavidPAClass &pa(*(DavidPAClass*)PairArray(typeNum));
    int ncomps=Path.Species(speciesNum).NumParticles;
    double lambda=Path.Species(speciesNum).lambda;
    for (int i=0;i<pa.kVals.size();i++){
      double k=pa.kVals(i);
      if (fequals(0.0,k,1e-10))
        yk_zero(speciesNum)=pa.uk_long(i);
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
           //cerr<<"VLONG IS "<<Vlong_k(speciesNum, j)<<" "<<pa.uk_long(i)<<" "<<speciesNum<<endl;
        }
      }
    }
    cout<<PathData.Path.Communicator.MyProc()<<" RPA Built: species("<<speciesNum<<") lambda("<<lambda<<") N("<<ncomps<<") Vol("<<vol<<")"<<endl;
    //for (int i=0;i<uk.extent(1);i++)
    //  cerr<<"KVecs: "<<sqrt(blitz::dot(Path.kVecs(i),Path.kVecs(i)))<<" "<<uk(speciesNum, i)<<" "<<duk(speciesNum, i)<<endl;
  }
}


void DavidLongRangeClassYk::ReadYk()
{
  for (int pai=0;pai<PairArray.size();pai++){
    DavidPAClass &pa(*((DavidPAClass*)PairArray(pai)));
    assert(pa.LongRangeDim==NDIM);
    for (int dim=0;dim<NDIM;dim++)
      assert(pa.LongRangeBox(dim)==Path.GetBox()[dim]);

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
    cout<<PathData.Path.Communicator.MyProc()<<" Masses: "<<pa.LongRangeMass1<<" "<<pa.LongRangeMass2<<" "<<Path.Species(specNum1).lambda<<" "<<Path.Species(specNum2).lambda<<endl;
    assert(fabs(pa.LongRangeMass1-Path.Species(specNum1).lambda)<1e-10);
    assert(fabs(pa.LongRangeMass2-Path.Species(specNum2).lambda)<1e-10);
  }
  //uk.resize(PairArray.size(), Path.kVecs.size());
  //duk.resize(PairArray.size(), Path.kVecs.size());
  //Vlong_k.resize(PairArray.size(), Path.kVecs.size());
  //yk_zero.resize(PairArray.size());
  uk.resize(Path.NumSpecies(),Path.kVecs.size());
  duk.resize(Path.NumSpecies(),Path.kVecs.size());
  Vlong_k.resize(Path.NumSpecies(),Path.kVecs.size());
  yk_zero.resize(Path.NumSpecies());
  if(PairArray.size()==1){
    cout<<PathData.Path.Communicator.MyProc()<<" Building RPA for single type" << endl;
    BuildRPA_SingleType();
  }
  else {
    cout <<PathData.Path.Communicator.MyProc()<< "Building RPA for multiple types" << endl;
    Build_MultipleSpecies();

  }
}

void DavidLongRangeClassYk::Read(IOSectionClass &in)
{
  assert(1==2);
  double myNum;
  uk.resize(Path.MagK.size());
  duk.resize(Path.MagK.size());
  for (int counter=0;counter<duk.size();counter++){
    duk(counter)=0.0;
  }
  string fileName;
  assert(in.ReadVar("LongRangeFile",fileName));
  ifstream infile;
  ///BUG: Currently hardcoded for actual file
  infile.open(fileName.c_str());
  string isRank;
  infile >> isRank;
  assert(isRank=="RANK");
  int is5; infile >> is5; assert(is5==5);
  int numkVec; infile >> numkVec;
  int is1; infile >>is1; assert(is1==1); infile >>is1; assert(is1==1);
  int is3; infile >> is3; assert(is3==3);
  int numLvl; infile >> numLvl; 
  infile >> isRank;
  assert(isRank=="BEGIN");
  infile >> isRank;
  assert(isRank=="k-space");
  infile >>isRank;
  assert(isRank=="action");

  cout<<"Loading David Long Range: "<<numLvl<<" "<<numkVec<<endl;
  for (int lvl=0;lvl<numLvl;lvl++)
    for (int isEnergy=0;isEnergy<3;isEnergy++)
      for (int kVec=0;kVec<numkVec;kVec++){
        infile>>myNum;
        if (lvl==0 && isEnergy==1){
          uk(kVec)=myNum;
        }
        if ((lvl==0 && isEnergy==2) || (lvl==0 && isEnergy==0)) {
          duk(kVec)+=myNum;
          cout<<"My energy is "<<myNum<<endl;
        }

      }

  infile.close();
}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}

inline double mag2 (const complex<double> &z1, const complex<double> &z2)
{
  return (z1.real()*z2.real() + z1.imag()*z2.imag());
}


/// Calculates the long range part of the action using David's breakup.
/// The short range part must be supplied as a dm file without the long
/// range part in it.  It ignores active particles.
double DavidLongRangeClassYk::SingleAction (int slice1, int slice2,const Array<int,1> &activeParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  int skip = 1<<level;

  //int species=0;
  set<int> speciesList;
  for(int p=0; p<activeParticles.size(); p++) {
    int ptcl = activeParticles(p);
    int spec = Path.ParticleSpeciesNum(ptcl);
    speciesList.insert(speciesList.begin(), spec);
  }

  bool only_do_inclusive=false;
  if (slice1==0 && slice2==PathData.Path.NumTimeSlices()-1)
    only_do_inclusive=false;
  else
    only_do_inclusive=true;
  //only_do_inclusive=false;
  if (GetMode() == NEWMODE){
    if (!only_do_inclusive)
      Path.UpdateRho_ks(slice1, slice2-skip, activeParticles, level);
    else
      Path.UpdateRho_ks(slice1+skip,slice2-skip,activeParticles,level);
  }

  double total=0;
  double factor;
  int startSlice=slice1;
  int endSlice=slice2-skip;
  if (only_do_inclusive){
    startSlice=slice1+skip;
    endSlice=slice2-skip;
  }

  //// Check to see if matrices are updated properly
  //Array<complex<double>,1> temp(Path.kVecs.size());
  //for (int slice=slice1; slice<=slice2; slice+=skip) {
  //  for (int species=0; species<Path.NumSpecies(); species++) {
  //  //       if (GetMode() == NEWMODE)
  //  //      Path.CalcRho_ks_Fast(slice,species);
  //  for (int ki=0; ki<Path.kVecs.size(); ki++)
  //    temp(ki) = Path.Rho_k(slice,species,ki);
  //    Path.CalcRho_ks_Fast(slice,species);
  //    for (int ki=0; ki<Path.kVecs.size(); ki++)
  //      if (mag2(temp(ki)-Path.Rho_k(slice,species,ki)) > 1.0e-12) 
  //      cerr << "Error in LongRangeClass::SingleAction.  "
  //      << "Cache inconsisency at slice=" 
  //      << slice << " species=" << Path.Species(species).Name 
  //      << "  Mode=" << GetMode() << endl
  //      << "temp  = " << temp(ki) << endl
  //      << "Rho_k = " << Path.Rho_k(slice,species,ki) << endl;
  //  }
  //}

  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=startSlice;slice<=endSlice;slice+=skip){
      // if ((slice == slice1) || (slice==slice2))
      //   factor = 0.5;
      // else
        factor = 1.0;
      for(set<int>::iterator it = speciesList.begin(); it!=speciesList.end(); it++) {
        int species = *it;
        double rhok2 = mag2(Path.Rho_k(slice,species,ki));
        total +=  factor*rhok2 * uk(species, ki);
      }
    }
  }

  // Cross-terms for Multiple Species
  for (int ki = 0; ki < Path.kVecs.size(); ki++) {
    for (int slice = startSlice; slice < endSlice; slice+=skip) {
      factor = 2.0;
      for (int species0 = 0; species0 < Path.NumSpecies()-1; species0++) {
        for (int species1 = species0+1; species1 < Path.NumSpecies(); species1++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          total += factor*rhok2 * uk(species0, ki);
        }
      }
    }
  }

  gettimeofday(&end,   &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return total;

}

/// Not really d_dbeta but total energy
double DavidLongRangeClassYk::d_dBeta (int slice1, int slice2,  int level)
{
  double total=0.0;
  double factor=1.0;
  for (int slice=slice1;slice<=slice2;slice++){
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
         double rhok2 = mag2(Path.Rho_k(slice,species,ki));
         sliceTotal +=  factor*rhok2 * duk(species, ki);
      }
    }
    total += sliceTotal;
  }

  // Cross-terms for Multiple Species
  for (int slice = slice1; slice <= slice2; slice++) {
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
          sliceTotal += factor*rhok2 * duk(species0, ki);
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
  for (int slice=slice1;slice<=slice2;slice++){
    double sliceTotal=0.0;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
         double rhok2 = mag2(Path.Rho_k(slice,species,ki));
         sliceTotal +=  factor*rhok2 * Vlong_k(species, ki);
      }
    }
    total += sliceTotal;
  }

  return total;

}



DavidLongRangeClassYk::DavidLongRangeClassYk(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix, Array<PairActionFitClass*,1> &pairArray)
  : ActionBaseClass (pathData), PairMatrix(pairMatrix), PairArray(pairArray)
{

}



string
DavidLongRangeClassYk::GetName()
{
  return "DavidsLongRange";
}
