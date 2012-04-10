#ifndef BSPLINE_HELPER_H
#define BSPLINE_HELPER_H

#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <blitz/tinyvec-et.h>
#include <complex>
#include "BoundaryCondition.h"
//#include "../Blitz.h"

using namespace blitz;
using namespace std;


// Solve for fixed first derivative.  
// data(0)  = dy_0/dx dx.  
// data(1)   = y_0 and
// data(M-1) = y_{M-1}.  
// data(M)   = dy_{M-1}/dx dx.  
template<typename T> inline void
SolveFirstDerivInterp1D (Array<T,1> &data, Array<T,1> &p)
{
  double ratio = 0.25;
  double ratio2 = 0.75;
  int M = data.size()-2;

  Array<double,1> d(M+2), mu(M+2);
  d = 1.5*data;
  mu = ratio;
  // First, eliminate leading coefficients
  mu(0) = 0.0;
  mu(M+1) = 0.0;
  mu(1) = ratio+ratio;
  d(0) /= (-ratio2);
  for (int row=1; row <=M; row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    mu(row) *= diagInv;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
  }
  d(M+1) /= -ratio2;
  mu(M+1) /= -ratio2;
  
  d(M+1)  -= d(M-1);
  mu(M+1) -= mu(M-1);
  d(M+1)  -= mu(M+1)*d(M);
  double diag = -1.0 - mu(M+1)*mu(M);
  p(M+1) = d(M+1)/diag;
 
  // Now go back upward, back substituting
  for (int row=M; row>=1; row--) 
    p(row) = d(row) - mu(row)*p(row+1);

  // And do 0th row
  p(0) = d(0) + p(2)/**d(2)*/;
}


// template<typename T> inline void
// SolvePeriodicInterp(Array<T,1> &data, Array<T,1> &p)
// {
//   double ratio = 0.25;
//   int N = data.size();

//   Array<double,1> d(N), gamma(N), mu(N);
//   d = 1.5*data;
//   p.resize(Range(0,N+2));
//   // First, eliminate leading coefficients
//   gamma (0) = ratio;
//   mu(0) = ratio;
//   mu(N-1) = ratio;
//   gamma(N-1) = 1.0;
//   for (int row=1; row <(N-1); row++) {
//     double diag = 1.0- mu(row-1)*ratio;
//     double diagInv = 1.0/diag;
//     gamma(row) = -ratio*gamma(row-1)*diagInv;
//     mu(row) = diagInv*ratio;
//     d(row)  = diagInv*(d(row)-ratio*d(row-1));
//     // Last row
//     d(N-1) -= mu(N-1) * d(row-1);
//     gamma(N-1) -= mu(N-1)*gamma(row-1);
//     mu(N-1) = -mu(N-1)*mu(row-1);
//   }
//   // Last row:  gamma(N-1) hold diagonal element
//   mu(N-1) += ratio;
//   gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
//   d(N-1) -= mu(N-1) * d(N-2);
//   p(N) = d(N-1)/gamma(N-1);
 
//   // Now go back upward, back substituting
//   for (int row=N-2; row>=0; row--) 
//     p(row+1) = d(row) - mu(row)*P(row+2) - gamma(row)*p(N);
// }



// 
template<typename T> inline void
SolvePeriodicInterp1D (Array<T,1> data, Array<T,1> p)
{
  double ratio = 0.25;
  int N = data.size();

  assert (p.size() == (N+3));

  Array<double,1> d(N), gamma(N), mu(N);
  d = 1.5*data;
  // First, eliminate leading coefficients
  gamma (0) = ratio;
  mu(0) = ratio;
  mu(N-1) = ratio;
  gamma(N-1) = 1.0;
  for (int row=1; row <(N-1); row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    gamma(row) = -ratio*gamma(row-1)*diagInv;
    mu(row) = diagInv*ratio;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
    // Last row
    d(N-1) -= mu(N-1) * d(row-1);
    gamma(N-1) -= mu(N-1)*gamma(row-1);
    mu(N-1) = -mu(N-1)*mu(row-1);
  }
  // Last row:  gamma(N-1) hold diagonal element
  mu(N-1) += ratio;
  gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
  d(N-1) -= mu(N-1) * d(N-2);
  p(N) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    p(row+1) = d(row) - mu(row)*p(row+2) - gamma(row)*p(N);
  p(0) = p(N);
  p(N+1) = p(1);
  p(N+2) = p(2);
}

inline void
SolvePeriodicInterp1D (Array<complex<double>,1> data, 
		       Array<complex<double>,1> p)
{
  int N = data.size();
  Array<double,1> dataReal(N), dataImag(N), pReal(N+3), pImag(N+3);
  for (int i=0; i<N; i++) {
    dataReal(i) = data(i).real();
    dataImag(i) = data(i).imag();
  }
  SolvePeriodicInterp1D(dataReal, pReal);
  SolvePeriodicInterp1D(dataImag, pImag);

  for (int i=0; i<N+3; i++) 
    p(i) = complex<double>(pReal(i), pImag(i));
}



// We multiply both sides by 1/4 to make diagonals equal 1
template<typename T> inline void
SolveDerivInterp1D (Array<T,1> data, Array<T,1> p,
		    TinyVector<double,4> abcdInitial,
		    TinyVector<double,4> abcdFinal)
{
  if (p.size() != (data.size()+2)) {
    cerr << "p.size() != (data.size()+2)\n";
    cerr << "p.size() = " << p.size() << endl;
    cerr << "data.size() = " << data.size() << endl;
    abort();
  }
  assert (p.size() == (data.size()+2));
  double al = 0.25*abcdInitial[0];
  double bl = 0.25*abcdInitial[1];
  double cl = 0.25*abcdInitial[2];
  double dl = 1.5 *abcdInitial[3];
  double ar = 0.25*abcdFinal[0];
  double br = 0.25*abcdFinal[1];
  double cr = 0.25*abcdFinal[2];
  double dr = 1.5 *abcdFinal[3];
    
  double ratio = 0.25;
  int M = data.size();

  Array<double,1> d(M+2), mu(M+2);
  d(Range(1,M)) = 1.5*data;
  mu = ratio;
  // First, eliminate leading coefficients
  double alInv = 1.0/al;
  bl *= alInv;
  cl *= alInv;
  dl *= alInv;

  d(0) = dl;  
  mu(0) = bl;
  mu(1) = ratio - ratio*cl;

  for (int row=1; row <=M; row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    mu(row) *= diagInv;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
  }
  
  br -= ar*mu(M-1);
  dr -= ar*d(M-1);
  cr -= br*mu(M);
  dr -= br*d(M);
  p(M+1) = dr/cr;
   
  // Now go back upward, back substituting
  for (int row=M; row>=1; row--) 
    p(row) = d(row) - mu(row)*p(row+1);

  // And do 0th row
  p(0) = dl -bl*p(1) - cl*p(2);
}


inline void
SolveDerivInterp1D (Array<complex<double>,1> data, 
		    Array<complex<double>,1> p,
		    TinyVector<double,4> abcdInitial,
		    TinyVector<double,4> abcdFinal)
{
  int N = data.size();
  int M = p.size();
  Array<double,1> dataReal(N), dataImag(N), pReal(M), pImag(M);
  for (int i=0; i<N; i++) {
    dataReal(i) = data(i).real();
    dataImag(i) = data(i).imag();
  }
  SolveDerivInterp1D(dataReal, pReal, abcdInitial, abcdFinal);
  SolveDerivInterp1D(dataImag, pImag, abcdInitial, abcdFinal);

  for (int i=0; i<M; i++) 
    p(i) = complex<double>(pReal(i), pImag(i));
}




#endif
