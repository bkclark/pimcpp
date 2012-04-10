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

//This is a comment designed to test the subversion submits. Feel free to remove it when you desire.
#include "PAclassicalFit.h"
#include "../Splines/BicubicSpline.h"

/// The following routines are used only if we are creating fits, not
/// using them.

void 
PAclassicalFitClass::ReadParams(IOSectionClass &inSection)
{
  //  UsePBC = inSection.ReadVar ("Box", Box);
}

void 
PAclassicalFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{ }


// void 
// PAclassicalFitClass::DoFit (Rho &rho)
// {
// }


// void 
// PAclassicalFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
// {
//   Uerror = 0.0;
//   dUerror = 0.0;
// }


void 
PAclassicalFitClass::WriteFit (IOSectionClass &outSection)
{
}



double 
PAclassicalFitClass::U(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  double V = Pot->V(q);
  //if (IsLongRange())
  //  V -= Vlong(q, level);
  return (beta*V);
}

double 
PAclassicalFitClass::dU(double q, double z, double s2, int level)
{
  double V = Pot->V(q);
  // if (IsLongRange())
  //  V -= Vlong(q, level);
  return (V);
}



bool 
PAclassicalFitClass::Read (IOSectionClass &in,
				double smallestBeta, int numBetas)
{
  SmallestBeta = smallestBeta;
  NumBetas=numBetas;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;
  assert (lambda == 0.0);

  // Read Potential;
  assert(in.OpenSection("Potential"));
  if (!in.ReadVar ("Z1Z2", Z1Z2))
    Z1Z2 = 0.0;
  Pot = ReadPotential(in);
  in.CloseSection();
  in.CloseSection();
  return true;
}



/////////////////////////
/// Long-ranged stuff ///
/////////////////////////
bool 
PAclassicalFitClass::IsLongRange()
{
  // This needs to be fixed.  We need to add this kind of function to
  // the potential base class. 
  return (Z1Z2 != 0.0);
}



double 
PAclassicalFitClass::Vlong(double q, int level)
{
  if (q <= 0.0)
    return 2.0/sqrt(M_PI)*Z1Z2*alpha;
  else 
    return Z1Z2/q*erf(alpha*q);
}

double 
PAclassicalFitClass::Vlong_k(double boxVol, double k, int level)
{
  if (k <= 0.0)
    k = 1.0e-30;
  double Vk =  4.0*M_PI*Z1Z2/(boxVol*k*k)*exp(-k*k/(4.0*alpha*alpha));
  //  cerr << "Vk = " << Vk << endl;
  return Vk;
}


// void PAclassicalFitClass::DoBreakup(const dVec& box,const Array<dVec,1> &kVecs)
// {
//     // Calculate the cutoff parameter
//   double minL, boxVol;
//   boxVol = minL = box[0];
//   for (int i=1; i<NDIM; i++) {
//     minL = min(minL, box[i]);
//     boxVol *= box[i];
//   }
//   alpha = 7.0/(minL*fabs(Z1Z2));
//   // Now, calculate the k-space parts
//   Ulong_k.resize(NumBetas, kVecs.size());
//   dUlong_k.resize(NumBetas, kVecs.size());
//   Ulong_0.resize(NumBetas);
//   dUlong_0.resize(NumBetas);
//   U_RPA_long_k.resize(NumBetas, kVecs.size());
//   dU_RPA_long_k.resize(NumBetas, kVecs.size());
//   for (int level=0; level<NumBetas; level++) {
//     double tau = SmallestBeta * pow(2.0, level);
//     Ulong_0(level)  = tau*Vlong(0.0, level);
//     dUlong_0(level) = Vlong(0.0, level);
//     for (int ki=0; ki<kVecs.size(); ki++) {
//       double k = sqrt(dot(kVecs(ki), kVecs(ki)));
//       Ulong_k(level,ki) = tau*Vlong_k(boxVol, k, level);
//       dUlong_k(level,ki) = Vlong_k(boxVol, k, level);
//     }
//   }
// }




/// The diagonal action only -- used for long-range breakup
double 
PAclassicalFitClass::Udiag(double q, int level)
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return beta*Pot->V(q);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::Udiag_p(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return beta*Pot->dVdr(q);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::Udiag_pp(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return beta*Pot->d2Vdr2(q);
}

/// The beta-derivative of the diagonal action
double 
PAclassicalFitClass::dUdiag    (double q, int level) 
{
  return Pot->V(q);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::dUdiag_p  (double q, int level) 
{
  return Pot->dVdr(q);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::dUdiag_pp (double q, int level) 
{
  return Pot->d2Vdr2(q);
}

/// The potential to which this action corresponds.
double 
PAclassicalFitClass::V  (double r) 
{
  return Pot->V(r);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::Vp (double r)
{
  return Pot->dVdr(r);
}

/// The q-derivative of the above
double 
PAclassicalFitClass::Vpp(double r) 
{
  return Pot->d2Vdr2(r);
}

/// HACK HACK HACK
/// This assumes we have a coulomb potential
double
PAclassicalFitClass::Xk_U(double k, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return -4.0*M_PI*beta*Z1Z2/(k*k) * cos(k*rCut);
}

double
PAclassicalFitClass::dXk_U_dk(double k, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return 4.0*M_PI*beta*Z1Z2/(k*k*k)*(2.0*cos(k*rCut)+k*rCut*sin(k*rCut));
}


/// HACK HACK HACK
/// This assumes we have a coulomb potential
double
PAclassicalFitClass::Xk_dU(double k, int level)
{
  return -4.0*M_PI*Z1Z2/(k*k) *cos(k*rCut);
}

/// HACK HACK HACK
/// This assumes we have a coulomb potential
double
PAclassicalFitClass::Xk_V(double k)
{
  return -4.0*M_PI*Z1Z2/(k*k) *cos(k*rCut);
}

/// HACK HACK HACK
/// This assumes we have a coulomb potential
double
PAclassicalFitClass::Vk (double k)
{
  return 4.0*M_PI*Z1Z2/(k*k);
}

void
PAclassicalFitClass::Setrc( double rc)
{
  rCut = rc;
}
