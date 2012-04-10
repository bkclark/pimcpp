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

#include <Common/MPI/Communication.h>
#include "LongRangeClass.h"
#include "../PathDataClass.h"
// #include <Common/Ewald/OptimizedBreakup.h>
// #include <Common/Integration/GKIntegration.h>

class CoulombXkIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  double k;
  double beta;
  JobType Task;

  inline double Uintegrand(double r)
  {
    double U = PA.Udiag(r, Level);
    U -= beta * PA.Z1Z2/r;
    return r * sin(k*r)*U;
  }
  inline double dUintegrand(double r)
  {
    double dU = PA.dUdiag(r, Level);
    dU -= PA.Z1Z2/r;
    return r * sin(k*r)*dU;
  }
  inline double Vintegrand(double r)
  {
    double V = PA.V(r);
    V -= PA.Z1Z2/r;
    return r * sin(k*r)*V;
  }
  
public:

  inline double operator()(double r)
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
  }

  CoulombXkIntegrand (PairActionFitClass &pa, int level, double k_,
		      JobType task) :
    PA(pa), Level(level), k(k_), Task(task)
  { 
    beta = pa.SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
  }
};

class XkIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  double k;
  JobType Task;
  inline double Uintegrand(double r)
  {
    double U = PA.Udiag(r, Level);
    return r * sin(k*r)*U;
  }
  inline double dUintegrand(double r)
  {
    double dU = PA.dUdiag(r, Level);
    return r * sin(k*r)*dU;
  }
  inline double Vintegrand(double r)
  {
    double V = PA.V(r);
    return r * sin(k*r)*V;
  } 

public:
  inline double operator()(double r) 
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  XkIntegrand (PairActionFitClass &pa, int level, double k_,
	       JobType task) :
    PA(pa), Level(level), k(k_), Task(task)
  { /* do nothing else*/  }
};



/// This calculates the quantity 
/// \f$ X_k \equiv -\frac{4 \pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
double LongRangeClass::CalcXk (int paIndex, int level, double k, double rc,
			    JobType task)
{
//   PathClass &Path = PathData.Path;
//   double absTol = 1.0e-7;
//   double relTol = 1.0e-5;

//   PairActionFitClass &pa = *PairArray(paIndex);
//   if (task == JOB_U)
//     absTol *= pa.SmallestBeta*pow(2.0, (double)level);

//   double Xk;
//   if (pa.Z1Z2 == 0.0) {
//     XkIntegrand integrand(pa, level, k, task);
//     GKIntegration<XkIntegrand, GK31> integrator(integrand);
//     Xk = 4.0*M_PI/(Path.GetVol()*k) * 
//       integrator.Integrate(rc, 20.0*rc, absTol, relTol, false);
//     return Xk;
//   }
//   else {
//     CoulombXkIntegrand integrand(pa, level, k, task);
//     GKIntegration<CoulombXkIntegrand, GK31> integrator(integrand);
//     // integrator.SetRelativeErrorMode();
    
    
//     if (false/*task != JOB_V*/) {
//       integrator.SetRelativeErrorMode();
//       Xk = -4.0*M_PI/(Path.GetVol()*k) * 
// 	integrator.Integrate (rc, 20.0*rc, relTol);    
//     }
//     else
//       Xk = -4.0*M_PI/(Path.GetVol()*k) * 
// 	integrator.Integrate(rc, 20.0*rc, absTol, relTol, false);
    
//     /// Add in the analytic part that I ignored
//     /// Multiply analytic term by tau only for U -- do not multiply
//     /// for dU or V.
    
//     double coef;
//     if (task == JOB_U) {
//       coef = pa.SmallestBeta;
//       for (int i=0; i<level; i++)
// 	coef *= 2.0;
//     }
//     else
//       coef = 1.0;

//     Xk -= coef*4.0*M_PI*pa.Z1Z2/(Path.GetVol()*k*k)*cos(k*rc);
//     return (Xk);
//   }
}


class UshortIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  JobType Task;
  inline double Uintegrand(double r)
  {
    double Ushort = PA.Udiag(r, Level) - PA.Ulong(Level)(r);
    return r*r*Ushort;
  }
  inline double dUintegrand(double r)
  {
    // This is a new attempt, using V instead of dU.
    double dUshort = PA.V(r) - PA.dUlong(Level)(r);
    return r*r*dUshort;
  }
  inline double Vintegrand(double r)
  {
    double Vshort = PA.V(r) - PA.Vlong(r);
    return r*r*Vshort;
  } 

public:
  inline double operator()(double r) 
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  UshortIntegrand (PairActionFitClass &pa, int level,
	       JobType task) :
    PA(pa), Level(level), Task(task)
  { /* do nothing else*/  }
};



class UlongIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  JobType Task;
  inline double Uintegrand(double r)
  {
    double Ulong = PA.Ulong(Level)(r);
    return r*r*Ulong;
  }
  inline double dUintegrand(double r)
  {
    double dUlong = PA.dUlong(Level)(r);
    return r*r*dUlong;
  }
  inline double Vintegrand(double r)
  {
    double Vlong = PA.Vlong(r);
    return r*r*Vlong;
  } 

public:
  inline double operator()(double r) 
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  UlongIntegrand (PairActionFitClass &pa, int level,
	       JobType task) :
    PA(pa), Level(level), Task(task)
  { /* do nothing else*/  }
};



void LongRangeClass::Read(IOSectionClass& in)
{
  //do nothing for now
}


LongRangeClass::LongRangeClass(PathDataClass &pathData,
			       Array<PairActionFitClass* ,2> &pairMatrix,
			       Array<PairActionFitClass*, 1> &pairArray) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix),
  PairArray(pairArray)
{
}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}

double 
LongRangeClass::SingleAction (int slice1, int slice2,
			      const Array<int,1> &changedParticles,
			      int level)
{
  if (GetMode() == NEWMODE)
    Path.UpdateRho_ks(slice1, slice2, changedParticles, level);
  
  int skip = (1<<level);
#ifdef DEBUG
  // Check to see if matrices are updated properly
  Array<complex<double>,1> temp(Path.kVecs.size());
  for (int slice=slice1; slice<=slice2; slice+=skip) {
    for (int species=0; species<Path.NumSpecies(); species++) {
//       if (GetMode() == NEWMODE)
// 	Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++)
	temp(ki) = Path.Rho_k(slice,species,ki);
      Path.CalcRho_ks_Fast(slice,species);
      for (int ki=0; ki<Path.kVecs.size(); ki++)
	if (mag2(temp(ki)-Path.Rho_k(slice,species,ki)) > 1.0e-12) 
	  cerr << "Error in LongRangeClass::SingleAction.  "
	       << "Cache inconsisency at slice=" 
	       << slice << " species=" << Path.Species(species).Name 
	       << "  Mode=" << GetMode() << endl
	       << "temp  = " << temp(ki) << endl
	       << "Rho_k = " << Path.Rho_k(slice,species,ki) << endl;
	  
      
    }
  }
#endif
      
  double homo = 0.0;
  double hetero = 0.0;
  double background = 0.0;
  double k0Homo = 0.0;
  double k0Hetero = 0.0;
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice+=skip) {
    cerr<<"doing slice "<<slice<<endl;
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    cerr<<"homologous "<<endl;
    for (int species=0; species<Path.NumSpecies(); species++) {
      // Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix(species,species);
      cerr<<"A"<<endl;
      if (pa.IsLongRange()) {
	cerr<<"B"<<endl;
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor * 0.5 * 2.0* rhok2 * pa.Ulong_k(level,ki);
	}
      }
      cerr<<"C"<<endl;
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      cerr<<"D"<<endl;
      homo -= factor * 0.5 * N * pa.Ulong_r0(level);
      // Or the neutralizing background term
      cerr<<"E"<<endl;
      background -= factor * 0.5*N*N*pa.Ushort_k0(level);
      cerr<<"F"<<endl;
      k0Homo += factor * 0.5*N*N*pa.Ulong_k0(level);
    }
    
    // Now do the heterologous terms
    cerr<<"heterologous "<<endl;
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
	    hetero += factor * 2.0 * rhorho * pa.Ulong_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  background  -= factor * N1*N2*pa.Ushort_k0(level);
	  k0Hetero += factor*N1*N2*pa.Ulong_k0(level);
	}
      }
    cerr<<"done slice "<<slice<<endl;
  }
  cerr<<"dout of loop"<<endl;
  double U = homo+hetero;
  if (UseBackground)
    U += background;
  else
    U += (k0Homo+k0Hetero);
  //  return (homo+hetero);
  return (U);
}


double LongRangeClass::d_dBeta (int slice1, int slice2,  int level)
{
  double homo = 0.0;
  double hetero = 0.0;
  double background = 0.0;
  double k0Homo = 0.0;
  double k0Hetero = 0.0;
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice+=skip) {
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &PA = *PairMatrix(species,species);
      if (PA.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor* 0.5 * 2.0* rhok2 * PA.dUlong_k(level,ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= factor * 0.5 * N * PA.dUlong_r0(level);
      // Or the neutralizing background term
      //background -= factor * 0.5*N*N*PA.dUshort_k0(level);
      // The background term is a constant and should be
      // independent of level.  Therefore, I'm using Vshort
      // instead of dUshort.
      background -= factor * 0.5*N*N*PA.dUshort_k0(level);
      // Or the k=0 terms
      k0Homo += factor*0.5*N*N*PA.dUlong_k0(level);
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &PA = *PairMatrix(species1,species2);
	if (PA.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += factor* 2.0 * rhorho * PA.dUlong_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  //	  background -= factor * N1*N2*PA.dUshort_k0(level);
	  // The background term is a constant and should be
	  // independent of level.  Therefore, I'm using Vshort
	  // instead of dUshort.
	  background -= factor * N1*N2*PA.dUshort_k0(level);
	  k0Hetero += factor*N1*N2*PA.dUlong_k0(level);
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

void
LongRangeClass::GradAction(int slice1, int slice2, 
			   const Array<int,1> &ptcls, int level,
			   Array<dVec,1> &gradVec)
{
  
  int skip = 1<<level;
  
  for (int slice = slice1; slice<=slice2; slice+=skip) {
    double factor = ((slice==slice1)||(slice==slice2)) ? 1.0 : 2.0;
    for (int pi=0; pi<ptcls.size(); pi++) {
      int ptcl = ptcls(pi);
      int species1 = Path.ParticleSpeciesNum(ptcl);
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	dVec r = Path(slice,ptcl);
	dVec k = Path.kVecs(ki);
	double phi = dot(r,k);
	double re, im;
	re=cos(phi);
	im=sin(phi);
	//	sincos(phi, &im, &re);
	complex<double> z(re, im);
	complex<double> rho_uSum(0.0, 0.0);
	for (int si=0; si<Path.NumSpecies(); si++) {
	  PairActionFitClass &PA = *PairMatrix(species1,si);
	  rho_uSum += PA.Ulong_k(level,ki) * Path.Rho_k(slice, si, ki);
	}
	rho_uSum = conj (rho_uSum);
	// Now, compute the imaginary part of the product, and
	// multiply by k
	gradVec(pi) -= 
	  factor*(z.real()*rho_uSum.imag() + z.imag()*rho_uSum.real())*k;
      }
    } // end pi loop
  } // end slice loop

}



void LongRangeClass::Init(IOSectionClass &in, IOSectionClass &out)
{
  if (PathData.Path.Getkc() == 0.0) {
    perr << "Missing kCutoff in System section.  Aborting.\n";
    abort();
  }
  int numKnots;
  assert(in.ReadVar ("NumBreakupKnots", numKnots));
  perr << "Doing optimized long range breakups...\n";
  if (PathData.Path.Communicator.MyProc() == 0) 
    out.NewSection ("LongRangeAction");
  OptimizedBreakup_U(numKnots, out);
  OptimizedBreakup_dU(numKnots, out);
  OptimizedBreakup_V(numKnots, out);
  if (PathData.Path.Communicator.MyProc() == 0) 
    out.CloseSection();
}




/// This computes the optimized breakups for the pair actions stored
/// in PairArray.  The parameters are the number of knots in
/// the "spline" representation of the long-range action and the
/// k-space cutoff.  
/// Only \f$\mathbf{k}\f$ with \f$|\mathbf{k}| < k_c$\f will be
/// included in the simulation sum.
void LongRangeClass::OptimizedBreakup_U(int numKnots, 
					IOSectionClass &out)
{
//   /// BUG: Long Range Optimized Breakups only work for NDIM=3
// #if NDIM==3
//   PathClass &Path = PathData.Path;
//   const double tolerance = 1.0e-7;
//   double kCut = Path.Getkc();
//   dVec box = Path.GetBox();
//   double boxVol = box[0]*box[1]*box[2];
//   double rc = 0.5*box[0];
//   for (int i=1; i<NDIM; i++)
//     rc = min (rc, 0.5*box[i]);
//   double kvol = Path.GetkBox()[0];
//   for (int i=1; i<NDIM; i++)
//     kvol *= Path.GetkBox()[i];
//   double kavg = pow(kvol,1.0/3.0);

  

//   LPQHI_BasisClass basis;
//   basis.Set_rc(rc);
//   basis.SetBox(box);
//   basis.SetNumKnots (numKnots);

//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double delta = basis.GetDelta();
//   double kMax = 20.0*M_PI/delta;
//   perr << "kCont = " << kCont 
//        << " kMax = " << kMax << endl;

//   OptimizedBreakupClass breakup(basis);
//   breakup.SetkVecs (kCut, kCont, kMax);
//   int numk = breakup.kpoints.size();
//   int N = basis.NumElements();
//   Array<double,1> t(N);
//   Array<bool,1>   adjust (N);
//   Array<double,1> Xk(numk);

//   // Would be 0.5, but with two timeslice distdisp, it could be a
//   // little longer
//   double rmax = 0.75 * sqrt (dot(box,box));
//   const int numPoints = 1000;
//   LongGrid.Init (0.0, rmax, numPoints);
//   Array<double,1> Ulong_r(numPoints), r(numPoints), Ushort_r(numPoints);
//   for (int i=0; i<numPoints; i++)
//     r(i) = LongGrid(i);

//   bool iAmRoot = PathData.Path.Communicator.MyProc() == 0;
//   if (iAmRoot) {
//     out.NewSection ("U");
//     out.WriteVar ("r", r);
//   }
//   for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
//     PairActionFitClass &pa = *PairArray(paIndex);
//     if (iAmRoot) {
//       perr << "Doing long range breakpus for species types (" 
// 	   << PairArray(paIndex)->Particle1.Name << ", " 
// 	   << PairArray(paIndex)->Particle2.Name << ")\n";
//       out.NewSection("PairAction");
//       out.WriteVar ("Particle1", pa.Particle1.Name);
//       out.WriteVar ("Particle2", pa.Particle2.Name);
//     }
//     pa.Setrc (rc);
//     pa.Ulong.resize(pa.NumBetas);
//     pa.Ulong_k.resize(pa.NumBetas,Path.kVecs.size());
//     pa.Ulong_k = 0.0;
//     pa.dUlong_dk.resize(pa.NumBetas,Path.kVecs.size());
//     pa.dUlong_dk = 0.0;
//     pa.Ulong_r0.resize(pa.NumBetas);
//     pa.Ushort_k0.resize(pa.NumBetas);
//     pa.Ulong_k0.resize(pa.NumBetas);
//     for (int level=0; level<pa.NumBetas; level++) {
//       if (iAmRoot)
// 	out.NewSection ("Level");
//       Ulong_r = 0.0;
      
//       // Calculate Xk's
//       perr << "Calculating Xk's for U...\n";
//       for (int ki=0; ki<numk; ki++) {
// 	//Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_U);
// 	//double oldXk = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_U);
// 	double k = breakup.kpoints(ki)[0];
// 	Xk(ki) = pa.Xk_U (k, level) / boxVol;
//       }
//       perr << "Done.\n";
      

//       // Set boundary conditions at rc:  Force value and first and
//       // second derivatives of long-range potential to match the full
//       // potential at rc.
//       adjust = true;
//       /// Warning:  the following constraints may cause instabilities!!!!
//       //t(N-3) = pa.Udiag(rc, level);                 adjust(N-3) = false;
//       //t(N-2) = pa.Udiag_p(rc, level)*delta;         adjust(N-2) = false;
//       //t(N-1) = pa.Udiag_pp(rc, level)*delta*delta;  adjust(N-1) = false;
//       //t(1) = 0.0;                                   adjust(1)   = false;
      
//       // Now, do the optimal breakup:  this gives me the coefficents
//       // of the basis functions, h_n in the array t.
//       perr << "Doing U breakup...\n";
//       breakup.DoBreakup (Xk, t, adjust);
//       perr << "Done.\n";
      
//       // Now, we must put this information into the pair action
//       // object.  First do real space part
//       pa.Ulong_r0(level)=0.0;
//       for (int n=0; n<N; n++)
// 	pa.Ulong_r0(level) += t(n)*basis.h(n,0.0);
//       for (int i=0; i<LongGrid.NumPoints; i++) {
// 	double r = LongGrid(i);
// 	if (r <= rc) {
// 	  // Sum over basis functions
// 	  for (int n=0; n<N; n++) 
// 	    Ulong_r(i) += t(n) * basis.h(n, r);
// 	}
// 	else
// 	  Ulong_r(i) = pa.Udiag (r, level);
//       }
//       pa.Ulong(level).Init(&LongGrid, Ulong_r);

//       if (iAmRoot) {
// 	// Now write to outfile
// 	out.WriteVar ("Ulong", Ulong_r);
// 	for (int i=0; i<numPoints; i++)
// 	  Ushort_r(i) = pa.Udiag(LongGrid(i), level) - Ulong_r(i);
// 	out.WriteVar ("Ushort", Ushort_r);
//       }
      

//       // Calculate FT of Ushort at k=0
//       UshortIntegrand shortIntegrand(pa, level, JOB_U);
//       GKIntegration<UshortIntegrand, GK31> shortIntegrator(shortIntegrand);
//       shortIntegrator.SetRelativeErrorMode();
//       pa.Ushort_k0(level) = 4.0*M_PI/boxVol * 
// 	shortIntegrator.Integrate(0.0, rc, tolerance);
//       perr << "Ushort_k0(" << level << ") = " << pa.Ushort_k0(level) << endl;

//       // Calculate FT of Ulong at k=0
//       UlongIntegrand longIntegrand(pa, level, JOB_U);
//       GKIntegration<UlongIntegrand, GK31> longIntegrator(longIntegrand);
//       longIntegrator.SetRelativeErrorMode();
//       pa.Ulong_k0(level) = 4.0*M_PI/boxVol * 
// 	longIntegrator.Integrate(0.0, rmax, tolerance);
//       perr << "Ulong_k0(" << level << ") = " << pa.Ulong_k0(level) << endl;


//       // Now do k-space part
//       for (int ki=0; ki < Path.kVecs.size(); ki++) {
// 	const dVec &kv = Path.kVecs(ki);
// 	double k = sqrt (dot(kv,kv));
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++) {
// 	  pa.Ulong_k(level,ki)   += t(n) * basis.c(n,k);
// 	  pa.dUlong_dk(level,ki) += t(n) * basis.dc_dk(n,k);
// 	}
// 	// Now add on part from rc to infinity
// 	// pa.Ulong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_U);
// 	pa.Ulong_k(level,ki)   -= pa.Xk_U  (k, level)/boxVol;
// #ifdef DEBUG
//  	double dXk_dk = pa.dXk_U_dk(k, level)/boxVol;
//  	double dXk_dkFD = 
//  	  (pa.Xk_U(k+1.0e-6,level)-pa.Xk_U(k-1.0e-6,level))/(2.0e-6*boxVol);
// 	if (PathData.Path.Communicator.MyProc() == 0)
// 	  fprintf (stderr, "analytic = %1.7e   FD = %1.7e\n",
// 		   dXk_dk, dXk_dkFD);
// #endif
// 	pa.dUlong_dk(level,ki) -= pa.dXk_U_dk(k, level)/boxVol;
//       }
// //       // HACK HACK HACK HACK
// //       FILE *fout = fopen ("Vlongk.dat", "w");
// //       for (double k=0; k<50.0; k+=0.01) {
// // 	double U = 0.0;
// // 	// Sum over basis functions
// // 	for (int n=0; n<N; n++)
// // 	  U += t(n) * basis.c(n,k);
// // 	fprintf (fout, "%1.16e %1.16e ", k, U);
// // 	// Now add on part from rc to infinity
// // 	U -= CalcXk(paIndex, level, k, rc);
// // 	fprintf (fout, "%1.16e \n", U);
// //       }
// //       fclose (fout);
//       if (iAmRoot)
// 	out.CloseSection (); // "Level"
//     }
//     if (iAmRoot)
//       out.CloseSection (); // "PairAction"
//   }
//   if (iAmRoot) {
//     out.CloseSection (); // "U"
//     out.FlushFile();
//   }
// #endif
}

void LongRangeClass::OptimizedBreakup_dU(int numKnots,
					 IOSectionClass &out)
{
//   bool iAmRoot = PathData.Path.Communicator.MyProc() == 0;

//   ///BUG: Optimized Breakup only works when NDIM==3
// #if NDIM==3
//   PathClass &Path = PathData.Path;
//   const double tolerance = 1.0e-7;
//   double kCut = Path.Getkc();
//   dVec box = Path.GetBox();
//   double boxVol = box[0]*box[1]*box[2];
//   double rc = 0.5*box[0];
//   for (int i=1; i<NDIM; i++)
//     rc = min (rc, 0.5*box[i]);
//   double kvol = Path.GetkBox()[0];
//   for (int i=1; i<NDIM; i++)
//     kvol *= Path.GetkBox()[i];
//   double kavg = pow(kvol,1.0/3.0);

//   LPQHI_BasisClass basis;
//   basis.Set_rc(rc);
//   basis.SetBox(box);
//   basis.SetNumKnots (numKnots);

//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double delta = basis.GetDelta();
//   double kMax = 20.0*M_PI/delta;
// //   perr << "kCont = " << kCont 
// //        << " kMax = " << kMax << endl;

//   OptimizedBreakupClass breakup(basis);
//   breakup.SetkVecs (kCut, kCont, kMax);
//   int numk = breakup.kpoints.size();
//   int N = basis.NumElements();
//   Array<double,1> t(N);
//   Array<bool,1>   adjust (N);
//   Array<double,1> Xk(numk);

//   // Would be 0.5, but with two timeslice distdisp, it could be a
//   // little longer
//   double rmax = 0.75 * sqrt (dot(box,box));
//   const int numPoints = 1000;
//   LongGrid.Init (0.0, rmax, numPoints);
//   Array<double,1> dUlong_r(numPoints), r(numPoints), dUshort_r(numPoints);
//   for (int i=0; i<numPoints; i++)
//     r(i) = LongGrid(i);

//   if (iAmRoot) {
//     out.NewSection ("dU");
//     out.WriteVar ("r", r);
//   }
//   for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
//     PairActionFitClass &pa = *PairArray(paIndex);
//     if (iAmRoot) {
//       out.NewSection("PairAction");
//       out.WriteVar ("Particle1", pa.Particle1.Name);
//       out.WriteVar ("Particle2", pa.Particle2.Name);
//     }

//     pa.Setrc (rc);
//     pa.dUlong.resize(pa.NumBetas);
//     pa.dUlong_k.resize(pa.NumBetas,Path.kVecs.size());
//     pa.dUlong_k = 0.0;
//     pa.dUlong_r0.resize(pa.NumBetas);
//     pa.dUshort_k0.resize(pa.NumBetas);
//     pa.dUlong_k0.resize(pa.NumBetas);
//     for (int level=0; level<pa.NumBetas; level++) {
//       if (iAmRoot) 
// 	out.NewSection ("Level");
//       dUlong_r = 0.0;

//       // Calculate Xk's
//       perr << "Calculating Xk's for dU...\n";
//       for (int ki=0; ki<numk; ki++) {	
// 	// Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_DU);
// 	double k = breakup.kpoints(ki)[0];
// 	Xk(ki) = pa.Xk_dU (k, level) / boxVol;
//       }
//       perr << "Done.\n";
//       // Set boundary conditions at rc:  Force value and first and
//       // second derivatives of long-range potential to match the full
//       // potential at rc.
//       adjust = true;
//       double delta = basis.GetDelta();
//       /// Warning:  the following constraints may cause instabilities!!!!
// //       t(N-3) = pa.dUdiag(rc, level);                 adjust(N-3) = false;
// //       t(N-2) = pa.dUdiag_p(rc, level)*delta;         adjust(N-2) = false;
// //       t(N-1) = pa.dUdiag_pp(rc, level)*delta*delta;  adjust(N-1) = false;
// //       t(1) = 0.0;                                    adjust(1)   = false;

//       // Now, do the optimal breakup:  this gives me the coefficents
//       // of the basis functions, h_n in the array t.
//       perr << "Doing dU breakup...\n";
//       breakup.DoBreakup (Xk, t, adjust);
//       perr << "Done.\n";
      
//       // Now, we must put this information into the pair action
//       // object.  First do real space part
//       pa.dUlong_r0(level)=0.0;
//       for (int n=0; n<N; n++)
// 	pa.dUlong_r0(level) += t(n)*basis.h(n,0.0);
//       for (int i=0; i<LongGrid.NumPoints; i++) {
// 	double r = LongGrid(i);
// 	if (r <= rc) {
// 	  // Sum over basis functions
// 	  for (int n=0; n<N; n++) 
// 	    dUlong_r(i) += t(n) * basis.h(n, r);
// 	}
// 	else
// 	  dUlong_r(i) = pa.dUdiag (r, level);
//       }
//       pa.dUlong(level).Init(&LongGrid, dUlong_r);
//       if (iAmRoot) {
// 	// Now write to outfile
// 	out.WriteVar ("dUlong", dUlong_r);
// 	for (int i=0; i<numPoints; i++)
// 	  dUshort_r(i) = pa.dUdiag(LongGrid(i), level) - dUlong_r(i);
// 	out.WriteVar ("dUshort", dUshort_r);
//       }


//       // Calculate FT of dUshort at k=0
//       UshortIntegrand shortIntegrand(pa, level, JOB_DU);
//       GKIntegration<UshortIntegrand, GK31> shortIntegrator(shortIntegrand);
//       shortIntegrator.SetRelativeErrorMode();
//       pa.dUshort_k0(level) = 4.0*M_PI/boxVol * 
// 	shortIntegrator.Integrate(0.0, rmax, tolerance);
//       perr << "dUshort_k0(" << level << ") = " << pa.dUshort_k0(level) << endl;

//       // Calculate FT of dUlong at k=0
//       UlongIntegrand longIntegrand(pa, level, JOB_DU);
//       GKIntegration<UlongIntegrand, GK31> longIntegrator(longIntegrand);
//       longIntegrator.SetRelativeErrorMode();
//       pa.dUlong_k0(level) = 4.0*M_PI/boxVol * 
// 	longIntegrator.Integrate(0.0, rmax, tolerance);
//       perr << "dUlong_k0(" << level << ") = " << pa.dUlong_k0(level) << endl;

//       // Now do k-space part
//       for (int ki=0; ki < Path.kVecs.size(); ki++) {
// 	const dVec &kv = Path.kVecs(ki);
// 	double k = sqrt (dot(kv,kv));
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++)
// 	  pa.dUlong_k(level,ki) += t(n) * basis.c(n,k);
// 	// Now add on part from rc to infinity
// 	//pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_DU);
// 	pa.dUlong_k(level,ki) -= pa.Xk_dU(k, level) / boxVol;
//       }
//       if (iAmRoot) 
// 	out.CloseSection (); // "Level"
//     }
//     if (iAmRoot) 
//       out.CloseSection (); // "PairAction"
//   }
//   if (iAmRoot) {
//     out.CloseSection (); // "dU"
//     out.FlushFile();
//   }
// #endif
}



// void LongRangeClass::OptimizedBreakup_dU(int numKnots)
// {
//   double kCut = Path.Getkc();
//   dVec box = Path.GetBox();
//   double rc = 0.5*box[0];
//   for (int i=1; i<NDIM; i++)
//     rc = min (rc, 0.5*box[i]);
//   double kvol = Path.GetkBox()[0];
//   for (int i=1; i<NDIM; i++)
//     kvol *= Path.GetkBox()[i];
//   double kavg = pow(kvol,1.0/3.0);
//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double kMax = 100 * kavg;
//   perr << "kCont = " << kCont 
//        << " kMax = " << kMax << endl;

//   LPQHI_BasisClass basis;
//   basis.Set_rc(rc);
//   basis.SetBox(box);
//   basis.SetNumKnots (numKnots);

//   OptimizedBreakupClass breakup(basis);
//   breakup.SetkVecs (kCut, kCont, kMax);
//   int numk = breakup.kpoints.size();
//   int N = basis.NumElements();
//   Array<double,1> t(N);
//   Array<bool,1>   adjust (N);
//   Array<double,1> Xk(numk);

//   // Would be 0.5, but with two timeslice distdisp, it could be a
//   // little longer
//   double rmax = 0.75 * sqrt (dot(box,box));
//   const int numPoints = 1000;
//   LongGrid.Init (0.0, rmax, numPoints);
//   Array<double,1> dUlong_r(numPoints);

//   for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
//     PairActionFitClass &pa = *PairArray(paIndex);
//     pa.dUlong.resize(pa.NumBetas);
//     pa.dUlong_k.resize(pa.NumBetas,Path.kVecs.size());
//     pa.dUlong_k = 0.0;
//     pa.dUlong_r0.resize(pa.NumBetas);
//     for (int level=0; level<pa.NumBetas; level++) {
//       dUlong_r = 0.0;
//       // Calculate Xk's
//       for (int ki=0; ki<numk; ki++)
// 	Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc,
// 			JOB_DU);

//       // Set boundary conditions at rc:  Force value and first and
//       // second derivatives of long-range potential to match the full
//       // potential at rc.
//       adjust = true;
//       double delta = basis.GetDelta();
//       perr << "dUdiag(rc) = " << pa.dUdiag(rc, level) << endl;
//       perr << "V(rc) = " << pa.V(rc) << endl;
//       perr << "dUdiag_p(rc) = " << pa.dUdiag_p(rc, level) << endl;
//       perr << "Vp(rc) = " << pa.Vp(rc) << endl;
//       perr << "dUdiag_pp(rc) = " << pa.dUdiag_pp(rc, level) << endl;
//       perr << "Vpp(rc) = " << pa.Vpp(rc) << endl;

      /// Warning:  the following constraints may cause instabilities!!!!
//       // t(N-3) = pa.dUdiag(rc, level);                 adjust(N-3) = false;
//       //t(N-2) = pa.dUdiag_p(rc, level)*delta;         adjust(N-2) = false;
//       //t(N-1) = pa.Vpp(rc)*delta*delta;  adjust(N-1) = false;
//       //t(1) = 0.0;                                    adjust(1)   = false;

// //       t(N-3) = pa.dUdiag(rc, level);     adjust(N-3) = false;
// //       t(N-2) = pa.dUdiag_p(rc, level);   adjust(N-2) = false;
// //       t(N-1) = pa.dUdiag_pp(rc, level);  adjust(N-1) = false;
// //       t(1) = 0.0;                        adjust(1)   = false;

//       // Now, do the optimal breakup:  this gives me the coefficents
//       // of the basis functions, h_n in the array t.
//       breakup.DoBreakup (Xk, t, adjust);
      
//       // Now, we must put this information into the pair action
//       // object.  First do real space part
//       pa.dUlong_r0(level)=0.0;
//       for (int n=0; n<N; n++)
// 	pa.dUlong_r0(level) += t(n)*basis.h(n,0.0);
//       for (int i=0; i<LongGrid.NumPoints; i++) {
// 	double r = LongGrid(i);
// 	if (r <= rc) {
// 	  // Sum over basis functions
// 	  for (int n=0; n<N; n++) 
// 	    dUlong_r(i) += t(n) * basis.h(n, r);
// 	}
// 	else
// 	  dUlong_r(i) = pa.dUdiag (r, level);
//       }
//       pa.dUlong(level).Init(&LongGrid, dUlong_r);

//       // Now do k-space part
//       for (int ki=0; ki < Path.kVecs.size(); ki++) {
// 	const dVec &kv = Path.kVecs(ki);
// 	double k = sqrt (dot(kv,kv));
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++)
// 	  pa.dUlong_k(level,ki) += t(n) * basis.c(n,k);
// 	// Now add on part from rc to infinity
// 	pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_DU);
//       }
//     }
//   }
// }


void LongRangeClass::OptimizedBreakup_V(int numKnots,
					IOSectionClass &out)
{
//   bool iAmRoot = PathData.Path.Communicator.MyProc() == 0;

//   ///BUG: Optimized breakup only works when NDIM==3
// #if NDIM==3
//   const double tolerance = 1.0e-7;
//   double kCut = Path.Getkc();
//   dVec box = Path.GetBox();
//   double boxVol = box[0]*box[1]*box[2];
//   double rc = 0.5*box[0];
//   for (int i=1; i<NDIM; i++)
//     rc = min (rc, 0.5*box[i]);
//   double kvol = Path.GetkBox()[0];
//   for (int i=1; i<NDIM; i++)
//     kvol *= Path.GetkBox()[i];
//   double kavg = pow(kvol,1.0/3.0);
// //   // We try to pick kcont to keep reasonable number of k-vectors
// //   double kCont = 50.0 * kavg;
// //   double kMax = 100 * kavg;
// //   perr << "kCont = " << kCont 
// //        << " kMax = " << kMax << endl;

//   LPQHI_BasisClass basis;
//   basis.Set_rc(rc);
//   basis.SetBox(box);
//   basis.SetNumKnots (numKnots);

//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double delta = basis.GetDelta();
//   double kMax = 20.0*M_PI/delta;
// //   perr << "kCont = " << kCont 
// //        << " kMax = " << kMax << endl;


//   OptimizedBreakupClass breakup(basis);
//   breakup.SetkVecs (kCut, kCont, kMax);
//   int numk = breakup.kpoints.size();
//   int N = basis.NumElements();
//   Array<double,1> t(N);
//   Array<bool,1>   adjust (N);
//   Array<double,1> Xk(numk);

//   // Would be 0.5, but with two timeslice distdisp, it could be a
//   // little longer
//   double rmax = 0.75 * sqrt (dot(box,box));
//   const int numPoints = 1000;
//   LongGrid.Init (0.0, rmax, numPoints);
//   Array<double,1> Vlong_r(numPoints), r(numPoints), Vshort_r(numPoints);
//   for (int i=0; i<numPoints; i++)
//     r(i) = LongGrid(i);


//   if (iAmRoot) {
//     out.NewSection ("V");
//     out.WriteVar ("r", r);
//   }

//   for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
//     PairActionFitClass &pa = *PairArray(paIndex);
//     if (iAmRoot) {
//       out.NewSection("PairAction");
//       out.WriteVar ("Particle1", pa.Particle1.Name);
//       out.WriteVar ("Particle2", pa.Particle2.Name);
//     }

//     pa.Setrc (rc);
//     pa.Vlong_k.resize(Path.kVecs.size());
//     pa.Vlong_k = 0.0;
//     Vlong_r = 0.0;

//     // Calculate Xk's
//     for (int ki=0; ki<numk; ki++) {
//       // Xk(ki) = CalcXk(paIndex, 0, breakup.kpoints(ki)[0], rc, JOB_V);
//       double k = breakup.kpoints(ki)[0];
//       Xk(ki) = pa.Xk_V (k) / boxVol;
//     }
    
//     // Set boundary conditions at rc:  Force value and first and
//     // second derivatives of long-range potential to match the full
//     // potential at rc.
//     adjust = true;
//     double delta = basis.GetDelta();
// //     t(N-3) = pa.V(rc);                 adjust(N-3) = false;
// //     t(N-2) = pa.Vp(rc)*delta;          adjust(N-2) = false;
// //     t(N-1) = pa.Vpp(rc)*delta*delta;   adjust(N-1) = false;
// //     t(1) = 0.0;                        adjust(1)   = false;
    
//     // Now, do the optimal breakup:  this gives me the coefficents
//     // of the basis functions, h_n in the array t.
//     perr << "Doing V breakup...\n";
//     breakup.DoBreakup (Xk, t, adjust);
//     perr << "Done.\n";
    
//     // Now, we must put this information into the pair action
//     // object.  First do real space part
//     pa.Vlong_r0=0.0;
//     for (int n=0; n<N; n++)
//       pa.Vlong_r0 += t(n)*basis.h(n,0.0);
//     for (int i=0; i<LongGrid.NumPoints; i++) {
//       double r = LongGrid(i);
//       if (r <= rc) {
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++) 
// 	  Vlong_r(i) += t(n) * basis.h(n, r);
//       }
//       else
// 	Vlong_r(i) = pa.V (r);
//     }
//     pa.Vlong.Init(&LongGrid, Vlong_r);
//     if (iAmRoot) {
//       // Now write to outfile
//       out.WriteVar ("Vlong", Vlong_r);
//       for (int i=0; i<numPoints; i++)
// 	Vshort_r(i) = pa.V(r(i)) - Vlong_r(i);
//       out.WriteVar ("Vshort", Vshort_r);
//     }

//     // Calculate FT of Ushort at k=0
//     UshortIntegrand shortIntegrand(pa, 0, JOB_V);
//     GKIntegration<UshortIntegrand, GK31> shortIntegrator(shortIntegrand);
//     shortIntegrator.SetRelativeErrorMode();
//     pa.Vshort_k0 = 4.0*M_PI/boxVol * 
//       shortIntegrator.Integrate(1.0e-100, rc, tolerance);
//     perr << "Vshort_k0 = " << pa.Vshort_k0 << endl;

//     // Calculate FT of Vlong at k=0
//     UlongIntegrand longIntegrand(pa, 0, JOB_V);
//     GKIntegration<UlongIntegrand, GK31> longIntegrator(longIntegrand);
//     longIntegrator.SetRelativeErrorMode();
//     pa.Vlong_k0 = 4.0*M_PI/boxVol * 
//       longIntegrator.Integrate(1.0e-100, rmax, tolerance);
//     perr << "Vlong_k0 = " << pa.Vlong_k0 << endl;


//     // Now do k-space part
//     for (int ki=0; ki < Path.kVecs.size(); ki++) {
//       const dVec &kv = Path.kVecs(ki);
//       double k = sqrt (dot(kv,kv));
//       // Sum over basis functions
//       for (int n=0; n<N; n++)
// 	pa.Vlong_k(ki) += t(n) * basis.c(n,k);
//       // Now add on part from rc to infinity
//       //pa.Vlong_k(ki) -= CalcXk(paIndex, 0, k, rc, JOB_V);
//       pa.Vlong_k(ki) -= pa.Xk_V(k) / boxVol;
//     }
//     if (iAmRoot)
//       out.CloseSection(); // "PairAction"
//   }

//   if (iAmRoot) {
//     out.CloseSection(); // V
//     out.FlushFile();
//   }
// #endif
}


string
LongRangeClass::GetName()
{
  return "LongRange";
}
