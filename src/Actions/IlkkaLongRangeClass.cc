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


#include "IlkkaLongRangeClass.h"
#include "../PathDataClass.h"
#include <set>
#include "sys/time.h"


bool IlkkaLongRangeClass::fequals(double a,double b,double tol)
{
  return abs(a-b)<tol;
}

bool IlkkaLongRangeClass::vecEquals(dVec &a, dVec &b,double tol)
{
  bool equals=true;
  for (int dim=0;dim<NDIM;dim++)
    equals = equals && fequals(a[dim],b[dim],tol);
  return equals;
}

void IlkkaLongRangeClass::Build_MultipleSpecies()
{
  uk.resize(PairArray.size(), Path.kVecs.size());
  duk.resize(PairArray.size(), Path.kVecs.size());
  Vlong_k.resize(PairArray.size(), Path.kVecs.size());
  uk0.resize(PairArray.size());
  duk0.resize(PairArray.size());
  vk0.resize(PairArray.size());
  ur0.resize(PairArray.size());
  dur0.resize(PairArray.size());
  vr0.resize(PairArray.size());

  double vol = Path.GetVol();
  for (int iPair=0; iPair<PairArray.size(); iPair++) {
    IlkkaPAClass &pa(*(IlkkaPAClass*)PairArray(iPair));

    // u
    for (int i=0; i<pa.k_u.size(); i++) {
      double k = pa.k_u(i);
      bool found = false;
      for (int j=0; j<Path.kVecs.size(); j++) {
        if (fequals(sqrt(dot(Path.kVecs(j),Path.kVecs(j))),k,1e-8)) {
          found = true;
          uk(iPair,j) = pa.uLong_k(i);
        }
      }
    }

    // du
    for (int i=0; i<pa.k_du.size(); i++) {
      double k = pa.k_du(i);
      bool found = false;
      for (int j=0; j<Path.kVecs.size(); j++) {
        if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-8)) {
          found = true;
          duk(iPair,j) = pa.duLong_k(i);
        }
      }
    }

    // v
    for (int i=0; i<pa.k_v.size(); i++) {
      double k = pa.k_v(i);
      bool found = false;
      for (int j=0; j<Path.kVecs.size(); j++) {
        if (fequals(sqrt(blitz::dot(Path.kVecs(j),Path.kVecs(j))),k,1e-8)) {
          found = true;
          Vlong_k(iPair,j) = pa.vLong_k(i);
        }
      }
    }


    // Constant terms
    uk0(iPair) = pa.uLong_k0;
    duk0(iPair) = pa.duLong_k0;
    ur0(iPair) = pa.uLong_r0;
    dur0(iPair) = pa.duLong_r0;
    vk0(iPair) = pa.vLong_k0;
    vr0(iPair) = pa.vLong_r0;

  }
}

void IlkkaLongRangeClass::WriteInfo(IOSectionClass &out)
{
  out.WriteVar ("Type", "Ilkka_Long_Range");
}

void IlkkaLongRangeClass::Read(IOSectionClass &in)
{
  TimeSpent=0;

  // Set flags
  Path.LongRange = true;
  Path.IlkkaLongRange = true;

  // Setup k Vecs and RhoK
  Path.SetupkVecs(in);

  // Get species number
  specNum1.resize(PairArray.size());
  specNum2.resize(PairArray.size());
  for (int iPair=0; iPair<PairArray.size(); iPair++){
    IlkkaPAClass &pa(*((IlkkaPAClass*)PairArray(iPair)));
    specNum1(iPair) = 0;
    while (Path.Species(specNum1(iPair)).Type!=pa.Particle1.Name){
      specNum1(iPair) += 1;
      assert(specNum1(iPair) < Path.NumSpecies());
    }
    specNum2(iPair) = 0;
    while (Path.Species(specNum2(iPair)).Type!=pa.Particle2.Name){
      specNum2(iPair) += 1;
      assert(specNum2(iPair) < Path.NumSpecies());
    }
  }

  // Build pair action arrays
  Build_MultipleSpecies();
}

/// Calculates the long range part of the action using Ilkka's breakup.
/// The short range part must be supplied as a dm file without the long
/// range part in it.  It ignores active particles.
double IlkkaLongRangeClass::SingleAction (int slice1, int slice2, const Array<int,1> &activeParticles, int level)
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

  double total = 0.;

  // Slower way
  //Path.DoPtcl = true;
  //for (int ptcl1Index=0; ptcl1Index<activeParticles.size(); ptcl1Index++) {
  //  int ptcl1 = activeParticles(ptcl1Index);
  //  Path.DoPtcl(ptcl1) = false;
  //  int species1 = Path.ParticleSpeciesNum(ptcl1);
  //  for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
  //    if (Path.DoPtcl(ptcl2)) {
  //      int species2 = Path.ParticleSpeciesNum(ptcl2);
  //      for (int ki=0; ki<Path.kVecs.size(); ki++) {
  //        for (int slice=startSlice; slice<=endSlice; slice++) {
  //          double factor = 2.0;
  //          dVec r;
  //          double rmag;
  //          PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
  //          double kr = dot(Path.kVecs(ki),r);
  //          total += factor * uk(PairIndex(species1,species2),ki) * cos(kr);
  //        }
  //      }
  //    }
  //  }
  //}

  // Homogolous terms
  for (int slice=startSlice; slice<=endSlice; slice+=skip) {
    for(set<int>::iterator it = speciesList.begin(); it!=speciesList.end(); it++) {
      int species = *it;
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
        double rhok2 = mag2(Path.Rho_k(slice,species,ki));
        total += 1.0*rhok2*uk(PairIndex(species,species),ki);
      }
    }
  }

  // Heterologous terms
  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=startSlice; slice<=endSlice; slice+=skip) {
      for (int species0=0; species0<Path.NumSpecies()-1; species0++) {
        for (int species1=species0+1; species1<Path.NumSpecies(); species1++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          total += 2.0*rhok2*uk(PairIndex(species0,species1),ki);
        }
      }
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return total;
}

double IlkkaLongRangeClass::d_dBeta (int slice1, int slice2,  int level)
{
  double total = 0.;
  double factor = 1.;

  // Slower way
  //Path.DoPtcl = true;
  //for (int ptcl1=0; ptcl1<Path.NumParticles(); ptcl1++) {
  //  Path.DoPtcl(ptcl1) = false;
  //  int species1 = Path.ParticleSpeciesNum(ptcl1);
  //  for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
  //    if (Path.DoPtcl(ptcl2)) {
  //      int species2 = Path.ParticleSpeciesNum(ptcl2);
  //      for (int ki=0; ki<Path.kVecs.size(); ki++) {
  //        for (int slice=slice1; slice<=slice2; slice++) {
  //          if ((slice==slice1) || (slice==slice2))
  //            factor = 1.0;
  //          else
  //            factor = 2.0;
  //          dVec r;
  //          double rmag;
  //          PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
  //          double kr = dot(Path.kVecs(ki),r);
  //          total += factor * duk(PairIndex(species1,species2),ki) * cos(kr);
  //        }
  //      }
  //    }
  //  }
  //}

  // Homogolous terms
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
         sliceTotal += factor*rhok2*duk(PairIndex(species,species),ki);
      }
    }
    total += sliceTotal;
  }

  //// Heterologous terms
  for (int slice=slice1; slice<=slice2; slice++) {
    double sliceTotal=0.0;
    if ((slice==slice1) || (slice==slice2))
      factor = 1.0;
    else
      factor = 2.0;
    for (int species0=0; species0<Path.NumSpecies()-1; species0++) {
      for (int species1=species0+1; species1<Path.NumSpecies(); species1++) {
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
double IlkkaLongRangeClass::V (int slice1, int slice2,  int level)
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
      for (int species1 = species0+1; species1 < Path.NumSpecies(); species1++) {
        for (int ki = 0; ki < Path.kVecs.size(); ki++) {
          double rhok2 = mag2(Path.Rho_k(slice,species0,ki),Path.Rho_k(slice,species1,ki));
          total += factor*rhok2*Vlong_k(PairIndex(species0,species1),ki);
        }
      }
    }
  }

  return total;

}

string IlkkaLongRangeClass::GetName()
{
  return "IlkkaLongRange";
}
