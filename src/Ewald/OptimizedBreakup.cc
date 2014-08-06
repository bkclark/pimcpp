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

// #include "../MPI/Communication.h"
#include "OptimizedBreakup.h"
#include "../MatrixOps/MatrixOps.h"
#include "../Integration/GKIntegration.h"

class cIntegrand
{
private:
  BasisClass &Basis;
  int n;
  double k;
public:
  double operator()(double r)
  {
    return 4.0*M_PI/(k*Basis.Omega) *Basis.h(n,r)*r*sin(k*r);
  }
  cIntegrand(BasisClass &basis, int n_, double k_) : Basis(basis), n(n_), k(k_)
  { /* do nothing */ }
};

///////////////////////////////////
/// BasisClass Member Functions ///
///////////////////////////////////
double BasisClass::c_numerical(int n, double k)
{
  int i = n/3;
  int numKnots = NumElements()/3;
  double ra = 0.0;
  double rb = r_c;

  cIntegrand integrand(*this, n, k);
  GKIntegration<cIntegrand,GK31> integrator(integrand);
  integrator.SetRelativeErrorMode();
  return integrator.Integrate (ra, rb, 1.0e-12, 1.0e-10, false);
}


void OptimizedBreakupClass::Addk(double k, double degeneracy)
{
  int ki=0;
  while ((ki < kpoints.size()) && (fabs(k-kpoints(ki)[0]) > 1.0e-12))
    ki++;
  if (ki == kpoints.size()) {
    kpoints.resizeAndPreserve(kpoints.size()+1);
    kpoints(ki)[0] = k;
    kpoints(ki)[1] = degeneracy;
  }
  else
    kpoints(ki)[1] += degeneracy;
}


///////////////////////////////////////////
/// OptimizedBreakupClass Memember Functions ///
///////////////////////////////////////////
void OptimizedBreakupClass::SetkVecs(double kc, double kCont, double kMax)
{
  double numk = 0;
  TinyVector<double,3> b;
  b[0] = 2.0*M_PI/Basis.GetBox()[0];
  b[1] = 2.0*M_PI/Basis.GetBox()[1];
  b[2] = 2.0*M_PI/Basis.GetBox()[2];
  TinyVector<int,3> maxIndex;
  maxIndex[0] = (int)ceil(kCont/b[0]);
  maxIndex[1] = (int)ceil(kCont/b[1]);
  maxIndex[2] = (int)ceil(kCont/b[2]);

  TinyVector<double,3> k;
  for (int ix=-maxIndex[0]; ix<=maxIndex[0]; ix++) {
    k[0] = ix*b[0];
    for (int iy=-maxIndex[1]; iy<=maxIndex[1]; iy++) {
      k[1] = iy*b[1];
      for (int iz=-maxIndex[2]; iz<=maxIndex[2]; iz++) {
	k[2] = iz*b[2];
	double k2 = dot(k,k);
	if ((k2 > (kc*kc)) && (k2 < (kCont*kCont))) {
	  Addk(sqrt(k2));
	  numk++;
	}
      }
    }
  }
  // Now, add kpoints to the list with approximate degeneracy.
  double kvol = b[0]*b[1]*b[2];
  const int N = 4000;
  double deltak = (kMax-kCont)/N;
  for (int i=0; i<N; i++) {
    double k1 = kCont+deltak*i;
    double k2 = k1+deltak;
    double k = 0.5*(k1+k2);
    double vol = 4.0*M_PI/3.0*(k2*k2*k2-k1*k1*k1);
    double degeneracy = vol/kvol;
    Addk(k, degeneracy);
    numk += degeneracy;
  }

  //cerr << "Total k vecs = " << numk << endl;
  //cerr << "non-degenerate k vecs = " << kpoints.size() << endl;
}

double OptimizedBreakupClass::DoBreakup(const Array<double,1> &Vk, 
					Array<double,1> &t)
{
  const double tolerance = 1.0e-16;
  //const double tolerance = 0.0;
  assert(t.rows()==Basis.NumElements());
  Array<double,2> A;
  Array<double,1> b;
  Array<double,2> cnk;

  int numElem = t.rows();
  A.resize(numElem, numElem);
  b.resize(numElem);
  cnk.resize(numElem,kpoints.rows());

  // Fill in cnk.
  for (int n=0; n<t.rows(); n++) {
    for (int ki=0; ki<kpoints.rows(); ki++) {
      double k = kpoints(ki)[0];
      cnk(n,ki) = Basis.c(n,k);
    }
  }

  // Now, fill in A and b
  A = 0.0;
  b = 0.0;
  for (int l=0; l<numElem; l++) {
    for (int ki=0; ki<kpoints.rows(); ki++) {
      b(l) += kpoints(ki)[1]*Vk(ki) * cnk(l, ki);
      for (int n=0; n<numElem; n++) 
	A(l,n) += kpoints(ki)[1]*cnk(l,ki)*cnk(n,ki);
    }
  }

  //  cerr << "A = " << A << endl;
  //cerr << "b = " << b << endl;

  // Now do SVD decomposition:
  Array<double,2> U(numElem, numElem), V(numElem, numElem);
  Array<double,1> S(numElem), Sinv(numElem);
  SVdecomp(A, U, S, V);
  
  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<S.size(); i++)
    Smax = max (S(i),Smax);

  for (int i=0; i<S.size(); i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv(i) == 0.0)
      numSingular++;
  if (numSingular > 0)
    cerr << "There were " << numSingular << " singular values.\n";
  t = 0.0;
  // Compute t_n, removing singular values
  for (int i=0; i<numElem; i++) {
    double coef = 0.0;
    for (int j=0; j<numElem; j++)
      coef += U(j,i) * b(j);
    coef *= Sinv(i);
    for (int k=0; k<numElem; k++)
      t(k) += coef * V(k,i);
  }
  // Calculate chi-squared
  double Yk, chi2;
  chi2 = 0.0;
  for (int ki=0; ki<kpoints.rows(); ki++) {
    Yk = Vk(ki);
    for (int n=0; n<t.rows(); n++) {
      Yk -= cnk(n,ki)*t(n);
    }
    chi2 += kpoints(ki)[1]*Yk+Yk;
  }
  return (chi2);
}

double OptimizedBreakupClass::DoBreakup(const Array<double,1> &Vk, Array<double,1> &t, const Array<bool,1> &adjust)
{
  const double tolerance = 1.0e-16;
  //const double tolerance = 0.0;
  assert(t.rows()==adjust.rows());
  assert(t.rows()==Basis.NumElements());
  Array<double,2> A;
  Array<double,1> b;
  Array<double,2> cnk;

  int N = t.rows();
  A.resize(N, N);
  b.resize(N);
  cnk.resize(N,kpoints.rows());

  // Fill in cnk.
  for (int n=0; n<t.rows(); n++) {
    for (int ki=0; ki<kpoints.rows(); ki++) {
      double k = kpoints(ki)[0];
      cnk(n,ki) = Basis.c(n,k);
    }
  }

  // Now, fill in A and b
  A = 0.0;
  b = 0.0;
  for (int l=0; l<N; l++) {
    for (int ki=0; ki<kpoints.rows(); ki++) {
      b(l) += kpoints(ki)[1]*Vk(ki) * cnk(l, ki);
      for (int n=0; n<N; n++) 
	A(l,n) += kpoints(ki)[1]*cnk(l,ki)*cnk(n,ki);
    }
  }
  // cerr << "A = " << A << endl;


  // Now reduce for constraints
  int M = N;
  for (int i=0; i<adjust.size(); i++) 
    if (!adjust(i))
      M--;

  // The c is for "constrained"
  Array<double,2> Ac(M,M);
  Array<double,1> bc(M), tc(M);

  // Build constrained Ac and bc
  int j=0;
  for (int col=0; col<N; col++) {
    if (adjust(col)) {
      // Copy column a A to Ac
      int i=0;
      for (int row=0; row<N; row++) 
	if (adjust(row)) {
	  Ac(i,j) = A(row,col);
	  i++;
	}
      j++;
    }
    else {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
      	b(row) -= A(row,col)*t(col);
    }
  }
  j=0;
  for (int row=0; row<N; row++)
    if (adjust(row)) {
      bc(j) = b(row);
      j++;
    }

  // Now do SVD decomposition:
  Array<double,2> U(M, M), V(M, M);
  Array<double,1> S(M), Sinv(M);
  SVdecomp(Ac, U, S, V);

  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<M; i++)
    Smax = max (S(i),Smax);

  for (int i=0; i<M; i++)
    if (S(i) < 0.0)
      cerr << "negative singlar value.\n";

  //  cerr << "Smax = " << Smax << endl;

  for (int i=0; i<M; i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));
  int numSingular = 0;
  for (int i=0; i<Sinv.size(); i++)
    if (Sinv(i) == 0.0)
      numSingular++;
  if (numSingular > 0)
    cerr << "There were " << numSingular << " singular values.\n";
  tc = 0.0;
  // Compute t_n, removing singular values
  for (int i=0; i<M; i++) {
    double coef = 0.0;
    for (int j=0; j<M; j++)
      coef += U(j,i) * bc(j);
    coef *= Sinv(i);
    for (int k=0; k<M; k++)
      tc(k) += coef * V(k,i);
  }

  // Now copy tc values into t
  j=0;
  for (int i=0; i<N; i++)
    if (adjust(i)) {
      t(i) = tc(j);
      j++;
    }

  // Calculate chi-squared
  double Yk, chi2;
  chi2 = 0.0;
  for (int ki=0; ki<kpoints.rows(); ki++) {
    Yk = Vk(ki);
    for (int n=0; n<t.rows(); n++) {
      Yk -= cnk(n,ki)*t(n);
    }
    chi2 += kpoints(ki)[1]*Yk*Yk;
  }
  return (chi2);
}

/////////////////////////////////////////
/// LPQHI_BasisClass Member Functions ///
/////////////////////////////////////////

void LPQHI_BasisClass::SetNumKnots(int n)
{
  assert (n > 1);
  NumKnots = n;
  if (r_c != 0.0) {
    delta = r_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
}

int LPQHI_BasisClass::NumElements()
{
  return 3*NumKnots;
}

void LPQHI_BasisClass::Set_rc(double rc)
{
  r_c = rc;
  if (NumKnots != 0) {
    delta = r_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
}


double LPQHI_BasisClass::h(int n, double r)
{
  int i=n/3;
  int alpha = n-3*i;
  double ra = delta*(i-1);
  double rb = delta*i;
  double rc = delta*(i+1);
  rc = min(r_c, rc);
  TinyVector<double,3> prefactor;
  prefactor =  1.0, 0.2, 0.02;
  if ((r > ra) && (r <= rb)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += (S(alpha,j) * prod);
      prod *= ((rb - r) * deltaInv);
    }
    //for (int j=0; j<alpha; j++)
    //  sum *= (-delta);
    for (int j=0; j<alpha; j++)
      sum *= -1.0;
    //sum *= prefactor[alpha];
    return (sum);
  }
  else if ((r > rb) && (r <= rc)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += S(alpha,j) * prod;
      prod *= ((r-rb) * deltaInv);
    }
    //for (int j=0; j<alpha; j++)
    //  sum *= delta;
    //sum *= prefactor[alpha];
    return sum;
  }
  return 0.0;
};


double LPQHI_BasisClass::c(int m, double k)
{
  int i=m/3;
  int alpha = m-3*i;
  
  double sum = 0.0;
  if (i == 0) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n));
    }
  else if (i == (NumKnots-1)) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dminus(i,k,n)*sign);
    }
  else
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n) + Dminus(i,k,n)*sign);
    }
  //for (int j=0; j<alpha; j++)
  //  sum *= delta;
  TinyVector<double,3> prefactor;
  prefactor = 1.0, 0.2, 0.02;
  //sum *= prefactor[alpha];
  return (sum);
};


double LPQHI_BasisClass::dc_dk(int m, double k)
{
  int i=m/3;
  int alpha = m-3*i;
  
  double sum = 0.0;
  if (i == 0) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (dDplus_dk(i,k,n));
    }
  else if (i == (NumKnots-1)) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (dDminus_dk(i,k,n)*sign);
    }
  else
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (dDplus_dk(i,k,n) + dDminus_dk(i,k,n)*sign);
    }
  //for (int j=0; j<alpha; j++)
  //  sum *= delta;
//   TinyVector<double,3> prefactor;
//   prefactor = 1.0, 0.2, 0.02;
  //sum *= prefactor[alpha];
  return (sum);
};
