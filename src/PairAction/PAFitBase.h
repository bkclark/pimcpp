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

#ifndef PA_FIT_BASE_H
#define PA_FIT_BASE_H

#include "../Communication/Communication.h" //include not needed
#include "../IO/IO.h"
#include "Potential.h"
#include "../Splines/CubicSpline.h"
#include "../Splines/BicubicSpline.h"
#include "../Splines/QuinticSpline.h"
#include "Particle.h"

/// Base class for pair actions
class PairActionFitClass
{
protected:
  void ReadHeader(IOSectionClass &inSection);
public:
  string Name;
  ParticleClass Particle1, Particle2;
  /// These store the long-range part of the action/potential in real
  /// space.  This will be subtracted from the total U to get the
  /// short-range part.
  blitz::Array<QuinticSpline,1> Ulong, dUlong;
  QuinticSpline Vlong;
  /// This stores the long-ranged part of the potential in k-space.  
  /// Indices: (level, k-point).  dUlong_k stores the beta-derivative.
  /// dUlong_dk stores the k-derivative of Ulong_k
  blitz::Array<double,2> Ulong_k, dUlong_k, dUlong_dk;
  blitz::Array<double,1> Vlong_k;
  /// This stores U_long(r=0);  Index is the level number.
  blitz::Array<double,1> Ulong_r0, dUlong_r0;
  double Vlong_r0;
  /// This stores (U/dU/V)_short(k=0);  Index is the level number.
  blitz::Array<double,1> Ushort_k0, dUshort_k0;
  double Vshort_k0;
  /// This stores (U/dU/V)_long(k=0);  Index is the level number.
  blitz::Array<double,1> Ulong_k0, dUlong_k0;
  double Vlong_k0;
  /// Stores the RPA form of the above.  This should be computed by
  /// ActionClass, since it couples all of the species pairs together
  /// and needs to know about the number of particles.
  blitz::Array<double,2> U_RPA_long_k, dU_RPA_long_k;
  /// This stores the beta-derivative of the above.
  bool SamplingTableRead;
  // Product of the two charges.  Zero if not coulomb or not charged.
  double Z1Z2;

  int NumBetas;
  
  Potential *Pot;
  double SmallestBeta;
  double lambda;

  virtual void ReadParams  (IOSectionClass &inSection) {};
  virtual void WriteBetaIndependentInfo (IOSectionClass &outSection) {};
  //  virtual void DoFit (Rho &rho) = 0;
  //  virtual void Error (Rho &rho, double &Uerror, double &dUerror)=0;
  virtual void WriteFit(IOSectionClass &outSection) {};

  virtual bool Read (IOSectionClass &inSection, double smallestBeta, int NumBetas) = 0;
  /// In the case of a long-ranged breakup, this should return Ushort
  virtual double U(double q, double z, double s2, int level) = 0;
  /// The beta-derivative of the action
  virtual double dU(double q, double z, double s2, int level) = 0;
  /// The potential to which this action corresponds.
  virtual double V  (double r) { return Pot->V(r); }
  /// The q-derivative of the above
  virtual double Vp (double r) { return Pot->dVdr(r); }
  /// The q-derivative of the above
  virtual double Vpp(double r) { return Pot->d2Vdr2(r); }
  /// These derivatives are needed to compute the gradient
  virtual void Derivs (double q, double z, double s2, int level,
		       double &d_dq, double &d_dz)
  { 
    cerr << "Error:  Derivs not implemented yet for this pair action type.";
    abort();
  }
  virtual void Derivs (double q, double z, double s2, int level,
		       double &d_dq, double &d_dz, double &d_ds)
  { 
    cerr << "Error:  Derivs not implemented yet for this pair action type.";
    abort();
  }

  /////////////////////////
  /// Long-ranged stuff ///
  /////////////////////////
  virtual bool IsLongRange() = 0;
  /// The diagonal action only -- used for long-range breakup
  virtual double Udiag(double q, int level)      { return 0.0; }
  /// The q-derivative of the above
  virtual double Udiag_p(double q, int level)    { return 0.0; }
  /// The q-derivative of the above
  virtual double Udiag_pp(double q, int level)   { return 0.0; }
  /// The beta-derivative of the diagonal action
  virtual double dUdiag    (double q, int level) { return 0.0; }
  /// The q-derivative of the above
  virtual double dUdiag_p  (double q, int level) { return 0.0; }
  /// The q-derivative of the above
  virtual double dUdiag_pp (double q, int level) { return 0.0; }
  /// This sets the cutoff radius for the long-range fit (NOT the core-radius for PH's)
  virtual void Setrc (double rc)                 {             }
  /** This is the k-compontent of the modified Fourier transform of
      the diagonal part of U, dU, and V, respectively, where the 
      integral is take from rc (not zero) to infinity. 
      \f[ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty \sin
      (kr) v(r) \, dr \f]
      Note that since we don't have the box here, we can't divide by
      \f$\Omega\f$. This must be done by the caller. */
  virtual double Xk_U  (double k, int level)     { return 0.0; }
  virtual double Xk_dU (double k, int level)     { return 0.0; }
  virtual double Xk_V  (double k)                { return 0.0; }
  virtual double Vk    (double k)                { return 0.0; }
  /// Returns the derivative w.r.t. k of Xk_U
  virtual double dXk_U_dk  (double k, int level) { return 0.0; }

  // Fills in the Vlong_k and dVlong_k array.
  virtual void DoBreakup (const Vec3 &box, const blitz::Array<Vec3,1> &kVecs) 
  { }
  PairActionFitClass() : Z1Z2(0.0), SamplingTableRead(false), Pot(NULL)
  { /* Do nothing */ }
};

#endif
