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
#include "LongRangePotClass.h"

LongRangePotClass::LongRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

inline double mag2(complex<double> z)
{
  return z.real()*z.real()+z.imag()*z.imag();
}

double LongRangePotClass::V(int slice)
{
  double homo = 0.0;
  double hetero = 0.0;
  double background = 0.0;
  double k0Homo = 0.0;
  double k0Hetero = 0.0;
  if (PathData.Actions.HaveLongRange()) {
    PathClass &Path = PathData.Path;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix(species,species);
      if (pa.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += 0.5 * 2.0* rhok2 * pa.Vlong_k(ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= 0.5 * N * pa.Vlong_r0;
      // Or the neutralizing background term
      background -= 0.5*N*N*pa.Vshort_k0;
      k0Homo += 0.5*N*N*pa.Vlong_k0;
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &pa= *PairMatrix(species1, species2);
	if (pa.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += 2.0 * rhorho * pa.Vlong_k(ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  background -= N1*N2*pa.Vshort_k0;
	  k0Hetero += N1*N2*pa.Vlong_k0;
	}
      }
  }
  double Vlong = homo+hetero;
  if (UseBackground)
    Vlong += background;
  else
    Vlong += (k0Homo+k0Hetero);
  return Vlong;
}
