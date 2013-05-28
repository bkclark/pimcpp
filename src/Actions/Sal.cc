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


#include "Sal.h"
#include "../PathDataClass.h"
#include <set>
#include  "sys/time.h"

void SalClass::Read(IOSectionClass &in)
{

}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}


/// Calculates the long range part of the action using David's breakup.
/// The short range part must be supplied as a dm file without the long
/// range part in it.  It ignores active particles.
double 
SalClass::SingleAction (int slice1, int slice2, 
			const Array<int,1> &activeParticles, 
			int level)

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

  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=startSlice;slice<=endSlice;slice+=skip){
      //      if ((slice == slice1) || (slice==slice2))
      //	factor = 0.5;
      //      else
	factor = 1.0;
      for(set<int>::iterator it = speciesList.begin(); it!=speciesList.end(); it++) {
        int species = *it;
        double rhok2 = mag2(Path.Rho_k(slice,species,ki));
        total +=  factor*rhok2;
      }
    }
  }

  gettimeofday(&end,   &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);


  return PathData.Path.tau*total;
    
}

  ///Not really d_dbeta but total energy
double SalClass::d_dBeta (int slice1, int slice2,  int level)
{
  double total=0.0;
  
  for (int ki=0; ki<Path.kVecs.size(); ki++) {
    for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++){
      double factor = 1.0;
      //      for(set<int>::iterator it = speciesList.begin(); it!=speciesList.end(); it++) {
      for (int species=0;species<PathData.Path.NumSpecies();species++){
        double rhok2 = mag2(Path.Rho_k(slice,species,ki));
        total +=  factor*rhok2;
      }
    }
  }
  return total;
}


  ///Not really d_dbeta but total energy
double SalClass::V (int slice1, int slice2,  int level)
{
  return 0.0;
}



SalClass::SalClass(PathDataClass &pathData):
  ActionBaseClass (pathData)
{

}



string
SalClass::GetName()
{
  return "Sal";
}
