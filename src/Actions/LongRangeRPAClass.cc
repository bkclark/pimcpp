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

#include "LongRangeRPAClass.h"

LongRangeRPAClass::LongRangeRPAClass(PathDataClass &pathData,
				     Array<PairActionFitClass* ,2> &pairMatrix,
				     Array<PairActionFitClass*, 1> &pairArray):
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix),
  PairArray(pairArray)
{

}

inline void order (int &species1, int &species2)
{
  if (species1 > species2) {
    int tmp = species1;
    species1 = species2;
    species2 = tmp;
  }
}


inline int uindex (int species1, int species2, int numSpecies)
{
  order(species1, species2);
  return (species2*(species2+1)/2 + species1);
}

inline int windex (int species1, int species2, int numSpecies)
{
  int offset = numSpecies*(numSpecies+1)/2;
  return (offset+uindex(species1, species2, numSpecies));
}


/// Note:  we must set Level and ki 
/// uwvec is indexed in the following way, 
/// given s1= species1 and s2 = species2, with s1<=s2
/// u_(s1,s2) = uwvec(
Array<double,1> LongRangeRPAClass::Integrand(double t, 
					     const Array<double,1> &uwvec)
{
  double levelTau = Path.tau;
  for (int i=0; i<Level; i++)
    levelTau *= 2.0;

  double boxVol = Path.GetVol();
  Array<double,1> duwvec(uwvec.size());
  int numSpecies = Path.NumSpecies();
  dVec &k = Path.kVecs(ki);
  double k2 = dot (k, k);
  // First, calculate the \f$\dot{u}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = 
	*PairMatrix(species1,species2);

      double vlongk;
      if (TaskIsU)
	vlongk = pa12.Ulong_k(Level, ki)/levelTau;
      else
	vlongk = pa12.dUlong_k(Level,ki);
      int i12 = uindex(species1,species2,numSpecies);
      double lambda12 = 
	0.5* (Path.Species(species1).lambda + Path.Species(species2).lambda);
      duwvec(i12) = -lambda12*k2*uwvec(i12) + vlongk;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double u23 = uwvec(uindex(species2, species3, numSpecies));
	double u13 = uwvec(uindex(species1, species3, numSpecies));
	double w23 = uwvec(windex(species2, species3, numSpecies));
	double w13 = uwvec(windex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*N3*lambda3*(u23*u13 + w23*u13 + u23*w13);
      }
    }
  // Next, calculate the \f$\dot{w}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = *PairMatrix(species1,species2);
      double k = sqrt(k2);
      double vlongk;
      if (TaskIsU)
	vlongk = pa12.Ulong_k(Level, ki)/levelTau;
      else
	vlongk = pa12.dUlong_k(Level,ki);
      double vshortk = pa12.Vk(k)/boxVol - vlongk;
      int i12 = windex(species1,species2,numSpecies);
      double lambda12 = 
	0.5* (Path.Species(species1).lambda + Path.Species(species2).lambda);
      duwvec(i12) = -lambda12*k2*uwvec(i12) + vshortk;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double w23 = uwvec(windex(species2, species3, numSpecies));
	double w13 = uwvec(windex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*N3*lambda3*(w23*w13);
      }
    }
  return duwvec;    
}

// All of the OptimizedBreakups must be computed before this is called.
void LongRangeRPAClass::Init(IOSectionClass &in)
{
  cerr << "Doing RPA correction...\n";

  for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
    PairActionFitClass &pa = *PairArray(paIndex);
    pa.U_RPA_long_k.resize(pa.NumBetas,Path.kVecs.size());
    pa.dU_RPA_long_k.resize(pa.NumBetas,Path.kVecs.size());
  }


  const int numPoints = 2000;
  double levelTau = Path.tau;
  int m = Path.NumSpecies()*(Path.NumSpecies()+1);
  Array<double,2> uwvec(numPoints, m);

  // Setup integrator
  //HACK  RungeKutta2<LongRangeRPAClass> integrator(*this);

  double boxVol = Path.GetVol();
  // Calculated RPA for U
  TaskIsU = true;
  ///For the RPA to work, every pairaction must have the same
  ///number of levels which is why the level is on the outside loop
  ///and we are using PairArray(0) as the thing that determines the
  ///number of levels
  for (Level=0; Level<PairArray(0)->NumBetas; Level++) {
    LinearGrid tGrid(0.0, levelTau, numPoints);
    for (ki=0; ki<Path.kVecs.size(); ki++) {
      //cerr << "ki = " << ki << endl;
      /// Set initial conditions
      for (int i=0; i<m; i++)
	uwvec(0, i) = 0.0;
      //HACK      integrator.Integrate(tGrid, 0, numPoints-1, uwvec);
      for (int species1=0; species1<Path.NumSpecies(); species1++)
	for (int species2=species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairMatrix(species1, species2);
	  pa.U_RPA_long_k(Level, ki) = 
	    uwvec(numPoints-1,uindex(species1, species2, Path.NumSpecies()));
	  double k2 = dot(Path.kVecs(ki), Path.kVecs(ki));
	  double k = sqrt(k2);
	  //cerr << "k = " << k << endl;
	  //cerr << "Numerical value = " << pa.U_RPA_long_k(Level,ki) << endl;
	  double lambda = pa.lambda;
	  double N = Path.Species(species1).NumParticles;
	  double vklong = pa.Ulong_k(Level, ki)/levelTau;
	  //cerr << "uklong  =         " << pa.Ulong_k(Level, ki) << endl;
	  //cerr << "vshort  =         " << 
	  //  pa.Vk(k)/boxVol - pa.Ulong_k(Level, ki)/levelTau << endl; 
	  //cerr << "ukshort =         " << 
	  //  uwvec(numPoints-1,windex(species1,species2,Path.NumSpecies()))
	  //     << endl;
	  double Uanalytic = (-1.0+sqrt(1.0+4.0*N*vklong/(lambda*k2)))/(2.0*N);
	  //cerr << "One component ground state analytic = " << Uanalytic 
	  //     << endl;
	}
    }
    levelTau *= 2.0;
  }

  levelTau = Path.tau;
  // Calculated RPA for dU
  TaskIsU = false;
  for (Level=0; Level<PairArray(0)->NumBetas; Level++) {
    LinearGrid tGrid(0.0, levelTau, numPoints);
    for (ki=0; ki<Path.kVecs.size(); ki++) {
      /// Set initial conditions
      for (int i=0; i<m; i++)
	uwvec(0, i) = 0.0;
      //HACK      integrator.Integrate(tGrid, 0, numPoints-1, uwvec);
      Array<double,1> duvec(m);
      duvec = Integrand (levelTau, uwvec(numPoints-1,Range::all()));
      for (int species1=0; species1<Path.NumSpecies(); species1++)
	for (int species2=species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairMatrix(species1, species2);
	  pa.dU_RPA_long_k(Level, ki) = 
	    duvec (uindex(species1, species2, Path.NumSpecies()));
	}
    } 
    levelTau *= 2.0;
  }
  Test();

  cerr << "done.\n";
  // Print out some debug info
  for (int i=0; i<PairArray.size(); i++) {
    PairActionFitClass& pa=*PairArray(i);
    string fname = pa.Particle1.Name + "-" +
      pa.Particle2.Name + ".dat";
    FILE *fout = fopen (fname.c_str(), "w");
    for (double q=0.0; q<20000.0; q+=1.0) {
      double Udiag = pa.Udiag(q, 0);
      double Ulong = pa.Ulong(0)(q);
      double dUdiag = pa.dUdiag(q, 0);
      double dUlong = pa.dUlong(0)(q);
      double V = pa.V(q);
      double Vlong = pa.Vlong(q);
      fprintf (fout, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", 
	       q, Udiag, Ulong, dUdiag, dUlong, V, Vlong);
    }
    fclose (fout);
  }
}


void LongRangeRPAClass::Test()
{
  FILE *fout = fopen ("RPAtest.dat", "w");
  double minL = Path.GetBox()[0];
  for (int i=1; i<NDIM; i++)
    minL = min(minL, Path.GetBox()[i]);
  int numPoints = 1000;
  LinearGrid rgrid(0.0, 0.5*minL, numPoints);
  int m = (Path.NumSpecies()*(Path.NumSpecies()+1))/2;
  Array<double,2> u(m,numPoints), du(m,numPoints), 
    uRPA(m,numPoints), duRPA(m, numPoints);
  u = 0.0; du = 0.0; uRPA = 0.0; duRPA = 0.0;
  for (int ri=0; ri<rgrid.NumPoints; ri++) {
    double r = rgrid(ri);
    dVec vr;
    vr[0] = r; vr[1] = 0.0; vr[2] = 0.0;
    for (int ki=0; ki<Path.kVecs.size(); ki++) {
      dVec &vk = Path.kVecs(ki);
      double coskr = cos(dot(vk, vr));
      
      for (int species1=0; species1<Path.NumSpecies(); species1++) 
	for (int species2 = species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairMatrix(species1, species2);
	  int si = uindex(species1, species2, Path.NumSpecies());
	  u(si, ri)     += pa.Ulong_k(0, ki) * 2.0*(1.0+coskr);
	  du(si, ri)    += pa.dUlong_k(0, ki) * 2.0*(1.0+coskr);
	  uRPA(si, ri)  += pa.U_RPA_long_k(0, ki) * 2.0*(1.0+coskr);
	  duRPA(si, ri) += pa.dU_RPA_long_k(0, ki) * 2.0*(1.0+coskr);
	}
    }
    fprintf (fout, "%1.12e", r);
    for (int si=0; si<m; si++) 
      fprintf (fout, " %1.12e %1.12e %1.12e %1.12e", 
	       u(si, ri), uRPA(si, ri), du(si, ri), duRPA(si,ri));
    fprintf (fout, "\n");
    
  }
  fclose (fout);
}


inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}
 

double 
LongRangeRPAClass::SingleAction (int slice1, int slice2,
				 const Array<int,1> &changedParticles,
				 int level)
{
  double homo = 0.0;
  double hetero = 0.0;
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice1+=skip) {
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix(species,species);
      if (pa.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor * 0.5 * 2.0* rhok2 * pa.U_RPA_long_k(level,ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= factor * 0.5 * N * pa.Ulong_r0(level);
      // Or the neutralizing background term
      //background -= factor * 0.5*N*N*pa.Ushort_k0(level);
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &pa = *PairMatrix(species1, species2);
	if (pa.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += factor * 2.0 * rhorho * pa.U_RPA_long_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  //background -= factor * N1*N2*pa.Ushort_k0(level);
	}
      }
  }
  double U;
  U = homo + hetero;

  return (homo+hetero);
}


double LongRangeRPAClass::d_dBeta (int slice1, int slice2,  int level)
{
  double homo = 0.0;
  double hetero = 0.0;
  double background = 0.0;
  double k0Homo = 0.0;
  double k0Hetero = 0.0;
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice1+=skip) {
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix (species, species);
      if (pa.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor* 0.5 * 2.0* rhok2 * pa.dU_RPA_long_k(level,ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= factor* 0.5 * N * pa.dUlong_r0(level);
      // Or the neutralizing background term
      background -= factor*0.5*N*N*pa.dUshort_k0(level);
      // Or the k=0 terms
      k0Homo += factor*0.5*N*N*pa.dUlong_k0(level);
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &pa = *PairMatrix(species1, species2);
	if (pa.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += factor* 2.0 * rhorho * pa.dU_RPA_long_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  background -= factor*N1*N2*pa.dUshort_k0(level);
	  k0Hetero += factor*N1*N2*pa.dUlong_k0(level);
	}
      }
  }

  double dU = homo+hetero;
  if (UseBackground)
    dU += background;
  else
    dU += (k0Homo+k0Hetero);
  return dU;
}

string
LongRangeRPAClass::GetName()
{
  return "LongRangeRPA";
}
