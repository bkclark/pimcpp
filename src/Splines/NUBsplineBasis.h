#ifndef NUB_SPLINE_BASIS_H
#define NUB_SPLINE_BASIS_H

#include "BoundaryCondition.h"
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <blitz/tinyvec-et.h>
#include <complex>
#include "Grid.h"
#ifdef __SSE3__
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

using namespace blitz;

template <typename GridType>
class NUBsplineBasis
{
private:
  GridType *GridPtr;

  // xVals is just the grid points, augmented by two extra points on
  // either side.  These are necessary to generate enough basis
  // functions. 
  Array<double,1> xVals;
  // dxInv(i)[j] = 1.0/(grid(i+j-1)-grid(i-2))
  // This is used to avoid division in evalating the splines
  Array<TinyVector<double,3>,1> dxInv;
  bool Periodic;
  
public:
  // Initialized xVals and dxInv.
  inline double Grid (int i) { return ((*GridPtr)(i)); }
  inline GridType& GetGrid() { return (*GridPtr); }
  inline void Init(GridType *gridPtr, bool periodic=false);
  // Evaluates the basis functions at a give value of x.  Returns the
  // index of the first basis function
  inline int  operator() (double x, TinyVector<double,4>& bfuncs) const;
  // Same as above, but also computes first derivatives
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv) const;
  // Same as above, but also computes second derivatives
  inline int  operator() (double x, TinyVector<float,4>& bfunc,
			  TinyVector<float,4> &deriv,
			  TinyVector<float,4> &deriv2) const;
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv,
			  TinyVector<double,4> &deriv2) const;
  // These versions take the grid point index.  These are needed for
  // the boundary conditions on the interpolating equations.  The
  // complex version return the basis functions in both the real and
  // imaginary parts.
  inline void operator() (int i, TinyVector<float,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs) const;

  inline void operator() (int i, TinyVector<float,4>& bfuncs,
			  TinyVector<float,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs) const;

  inline void operator() (int i, TinyVector<float,4>& bfuncs,
			  TinyVector<float,4> &dbfuncs,
			  TinyVector<float,4> &d2bfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs,
			  TinyVector<double,4> &d2bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs,
			  TinyVector<complex<double>,4> &d2bfuncs) const;
};

template<typename GridType> void
NUBsplineBasis<GridType>::Init(GridType *gridPtr, bool periodic)
{
  Periodic = periodic;
  GridPtr = gridPtr;
  GridType &grid =  (*GridPtr);

  int N = grid.NumPoints;
  dxInv.resize(N+3);
  xVals.resize(N+6);
  for (int i=0; i<N; i++)
    xVals(i+2) = grid(i);
  if (!Periodic) {
    xVals(0) = grid(0) - 2.0*(grid(1)-grid(0));
    xVals(1) = grid(0) - 1.0*(grid(1)-grid(0));
    xVals(N+2) = grid(N-1) + 1.0*(grid(N-1)-grid(N-2));
    xVals(N+3) = grid(N-1) + 2.0*(grid(N-1)-grid(N-2));
    xVals(N+4) = grid(N-1) + 3.0*(grid(N-1)-grid(N-2));
    xVals(N+5) = grid(N-1) + 4.0*(grid(N-1)-grid(N-2));
  }
  else {
    xVals(1)   = grid(0)   - (grid(N-1) - grid(N-2));
    xVals(0)   = grid(0)   - (grid(N-1) - grid(N-3));
    xVals(N+2) = grid(N-1) + (grid(1)   - grid(0));
    xVals(N+3) = grid(N-1) + (grid(2)   - grid(0));
    xVals(N+4) = grid(N-1) + (grid(3)   - grid(0));
  }
  for (int i=0; i<N+3; i++) 
    for (int j=0; j<3; j++) 
      dxInv(i)[j] = 1.0/(xVals(i+j+1)-xVals(i));
}


template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, TinyVector<double,4> &bfuncs) const
{
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  int i = grid.ReverseMap (x);
  int i2 = i+2;
//   b1[0] = (grid(i+1)-x)/(grid(i+1)-grid(i));
//   b1[1] = (x-grid(i))/(grid(i+1)-grid(i));
//   b2[0] = (grid(i+1)-x)/(grid(i+1)-grid(i-1))*b1[0];
//   b2[1] = ((x-grid(i-1))/(grid(i+1)-grid(i-1))*b1[0] +
// 	   (grid(i+2)-x)/(grid(i+2)-grid(i))*b1[1]);
//   b2[2] = (x-grid(i))/(grid(i+2)-grid(i))*b1[1];
//   bfuncs[0] = (grid(i+1)-x)/(grid(i+1)-grid(i-2)) * b2[0];
//   bfuncs[1] = ((x-grid(i-2))/(grid(i+1)-grid(i-2))*b2[0] +
// 	       (grid(i+2)-x)/(grid(i+2)-grid(i-1))*b2[1]);
//   bfuncs[2] = ((x-grid(i-1))/(grid(i+2)-grid(i-1))*b2[1] +
// 	       (grid(i+3)-x)/(grid(i+3)-grid(i))*b2[2]);
//   bfuncs[3] = (x-grid(i))/(grid(i+3)-grid(i))*b2[2];

//   b1[0]     = (grid(i+1)-x)  * dxInv(i+2)[0];
//   b1[1]     = (x-grid(i))    * dxInv(i+2)[0];
//   b2[0]     = (grid(i+1)-x)  * dxInv(i+1)[1] * b1[0];
//   b2[1]     = ((x-grid(i-1)) * dxInv(i+1)[1] * b1[0]+
// 	       (grid(i+2)-x) * dxInv(i+2)[1] * b1[1]);
//   b2[2]     = (x-grid(i))    * dxInv(i+2)[1] * b1[1];
//   bfuncs[0] = (grid(i+1)-x)  * dxInv(i  )[2] * b2[0];
//   bfuncs[1] = ((x-grid(i-2)) * dxInv(i  )[2] * b2[0] +
// 	       (grid(i+2)-x) * dxInv(i+1)[2] * b2[1]);
//   bfuncs[2] = ((x-grid(i-1)) * dxInv(i+1)[2] * b2[1] +
// 	       (grid(i+3)-x) * dxInv(i+2)[2] * b2[2]);
//   bfuncs[3] = (x-grid(i))    * dxInv(i+2)[2] * b2[2];

  b1[0]     = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]     = (x-xVals(i2))    * dxInv(i+2)[0];
  b2[0]     = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]     = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
	       (xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]     = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  bfuncs[0] = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
	       (xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
	       (xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3] = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];
  return i;
}


template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  double x = grid(i);

  b1[0]     = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]     = (x-xVals(i2))    * dxInv(i+2)[0];
  b2[0]     = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]     = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
	       (xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]     = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  bfuncs[0] = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
	       (xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
	       (xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3] = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i,
				     TinyVector<complex<double>,4> &bfuncs) const
{
  TinyVector<double,4> funcs;
  (*this)(i, funcs);
  bfuncs[0] = complex<double>(funcs[0], funcs[0]);
  bfuncs[1] = complex<double>(funcs[1], funcs[1]);
  bfuncs[2] = complex<double>(funcs[2], funcs[2]);
  bfuncs[3] = complex<double>(funcs[3], funcs[3]);
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i,
				     TinyVector<float,4> &bfuncs) const
{
  TinyVector<double,4> funcs;
  (*this)(i, funcs);
  bfuncs[0] = float(funcs[0]);
  bfuncs[1] = float(funcs[1]);
  bfuncs[2] = float(funcs[2]);
  bfuncs[3] = float(funcs[3]);
}


template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs) const
{
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  int i = grid.ReverseMap (x);
  int i2 = i+2;

  b1[0]      = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]      = (x-xVals(i2))    * dxInv(i+2)[0];
  
  b2[0]      = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]      = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
		(xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]      = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  
  bfuncs[0]  = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1]  = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
		(xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2]  = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
		(xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3]  = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  return i;


  // The following is the stupid way of doing things:
  //   double  db1[2], db2[3];
  //   db1[0]     = -dxInv(i+2)[0];
  //   db1[1]     =  dxInv(i+2)[0];
  //   db2[0]     = ((xVals(i2+1)-x)  * dxInv(i+1)[1] * db1[0] - dxInv(i+1)[1] * b1[0]);
  //   db2[1]     = (((x-xVals(i2-1)) * dxInv(i+1)[1] * db1[0] + dxInv(i+1)[1] * b1[0]) +
  // 		((xVals(i2+2)-x) * dxInv(i+2)[1] * db1[1] - dxInv(i+2)[1] * b1[1]));
  //   db2[2]     = ((x-xVals(i2))    * dxInv(i+2)[1] * db1[1] + dxInv(i+2)[1] * b1[1]);
  //   dbfuncs[0] = ((xVals(i2+1)-x) * dxInv(i  )[2] * db2[0] - dxInv(i  )[2] * b2[0]);
  //   dbfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * db2[0] + dxInv(i  )[2] * b2[0] +
  // 		(xVals(i2+2)-x) * dxInv(i+1)[2] * db2[1] - dxInv(i+1)[2] * b2[1]);
  //   dbfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * db2[1] + dxInv(i+1)[2] * b2[1] +
  // 		(xVals(i2+3)-x) * dxInv(i+2)[2] * db2[2] - dxInv(i+2)[2] * b2[2]);
  //   dbfuncs[3] = ((x-xVals(i2))   * dxInv(i+2)[2] * db2[2] + dxInv(i+2)[2] * b2[2]);
}



template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double x = grid(i);
  double b1[2], b2[3];

  b1[0]     = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]     = (x-xVals(i2))    * dxInv(i+2)[0];
  b2[0]     = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]     = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
	       (xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]     = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  bfuncs[0] = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
	       (xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
	       (xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3] = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				     TinyVector<complex<double>,4> &dbfuncs) const
{
  TinyVector<double,4> funcs, dfuncs;
  (*this)(i, funcs, dfuncs);
  bfuncs[0]  = complex<double>(funcs[0],  funcs[0]);
  bfuncs[1]  = complex<double>(funcs[1],  funcs[1]);
  bfuncs[2]  = complex<double>(funcs[2],  funcs[2]);
  bfuncs[3]  = complex<double>(funcs[3],  funcs[3]);
  dbfuncs[0] = complex<double>(dfuncs[0], dfuncs[0]);
  dbfuncs[1] = complex<double>(dfuncs[1], dfuncs[1]);
  dbfuncs[2] = complex<double>(dfuncs[2], dfuncs[2]);
  dbfuncs[3] = complex<double>(dfuncs[3], dfuncs[3]);
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<float,4> &bfuncs,
				     TinyVector<float,4> &dbfuncs) const
{
  TinyVector<double,4> funcs, dfuncs;
  (*this)(i, funcs, dfuncs);
  bfuncs[0]  = float(funcs[0] );
  bfuncs[1]  = float(funcs[1] );
  bfuncs[2]  = float(funcs[2] );
  bfuncs[3]  = float(funcs[3] );
  dbfuncs[0] = float(dfuncs[0]);
  dbfuncs[1] = float(dfuncs[1]);
  dbfuncs[2] = float(dfuncs[2]);
  dbfuncs[3] = float(dfuncs[3]);
}

template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, 
				     TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs,
				     TinyVector<double,4> &d2bfuncs) const
{
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  int i = min(grid.ReverseMap (x), grid.NumPoints-1);
  int i2 = i+2;

  b1[0]      = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]      = (x-xVals(i2))    * dxInv(i+2)[0];
  
  b2[0]      = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]      = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
		(xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]      = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  
  bfuncs[0]  = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1]  = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
		(xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2]  = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
		(xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3]  = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  d2bfuncs[0] = 6.0 * (+dxInv(i+0)[2]*dxInv(i+1)[1]*b1[0]);
  d2bfuncs[1] = 6.0 * (-dxInv(i+1)[1]*(dxInv(i+0)[2]+dxInv(i+1)[2])*b1[0] +
		       dxInv(i+1)[2]*dxInv(i+2)[1]*b1[1]);
  d2bfuncs[2] = 6.0 * (dxInv(i+1)[2]*dxInv(i+1)[1]*b1[0] -
		       dxInv(i+2)[1]*(dxInv(i+1)[2] + dxInv(i+2)[2])*b1[1]);
  d2bfuncs[3] = 6.0 * (+dxInv(i+2)[2]*dxInv(i+2)[1]*b1[1]);

  return i;
}



template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				     TinyVector<complex<double>,4> &dbfuncs,
				     TinyVector<complex<double>,4> &d2bfuncs) const
{
  TinyVector<double,4> funcs, dfuncs, d2funcs;
  (*this)(i, funcs, dfuncs, d2funcs);
  bfuncs[0]   = complex<double>(  funcs[0],   funcs[0]);
  bfuncs[1]   = complex<double>(  funcs[1],   funcs[1]);
  bfuncs[2]   = complex<double>(  funcs[2],   funcs[2]);
  bfuncs[3]   = complex<double>(  funcs[3],   funcs[3]);
  dbfuncs[0]  = complex<double>( dfuncs[0],  dfuncs[0]);
  dbfuncs[1]  = complex<double>( dfuncs[1],  dfuncs[1]);
  dbfuncs[2]  = complex<double>( dfuncs[2],  dfuncs[2]);
  dbfuncs[3]  = complex<double>( dfuncs[3],  dfuncs[3]);
  d2bfuncs[0] = complex<double>(d2funcs[0], d2funcs[0]);
  d2bfuncs[1] = complex<double>(d2funcs[1], d2funcs[1]);
  d2bfuncs[2] = complex<double>(d2funcs[2], d2funcs[2]);
  d2bfuncs[3] = complex<double>(d2funcs[3], d2funcs[3]);
}



template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs,
				     TinyVector<double,4> &d2bfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double x = grid(i);
  double b1[2], b2[3];

  b1[0]     = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]     = (x-xVals(i2))    * dxInv(i+2)[0];
  b2[0]     = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]     = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
	       (xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]     = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  bfuncs[0] = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
	       (xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
	       (xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3] = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  d2bfuncs[0] = 6.0 * (+dxInv(i+0)[2]*dxInv(i+1)[1]*b1[0]);
  d2bfuncs[1] = 6.0 * (-dxInv(i+1)[1]*(dxInv(i+0)[2]+dxInv(i+1)[2])*b1[0] +
		       dxInv(i+1)[2]*dxInv(i+2)[1]*b1[1]);
  d2bfuncs[2] = 6.0 * (dxInv(i+1)[2]*dxInv(i+1)[1]*b1[0] -
		       dxInv(i+2)[1]*(dxInv(i+1)[2] + dxInv(i+2)[2])*b1[1]);
  d2bfuncs[3] = 6.0 * (+dxInv(i+2)[2]*dxInv(i+2)[1]*b1[1]);
}

#ifdef __SSE3__
  #include <xmmintrin.h>
#endif

////////////////////////////////////////////////////////////
//              LinearGrid specialization                 //
////////////////////////////////////////////////////////////
template<>
class NUBsplineBasis<LinearGrid>
{
private:
  LinearGrid *GridPtr;
  TinyMatrix<double,4,4> A, dA, d2A, d3A;
  double GridStart, GridEnd, GridDelta, GridDeltaInv, L, Linv;
  double Periodic, Sixth, TwoThirds;
  inline int Find (double x) const;
  mutable TinyVector<double,4> tp;
#ifdef __SSE3__
  mutable __m128 _tp;
  mutable __m128 _A[4], _dA[4], _d2A[4], _GDE;
#endif
  mutable int i0;
public:
  // Initialized xVals and dxInv.
  inline double Grid (int i) { return ((*GridPtr)(i)); }
  inline LinearGrid& GetGrid() { return (*GridPtr); }
  inline void Init(LinearGrid *gridPtr, bool periodic=false);
  // Evaluates the basis functions at a give value of x.  Returns the
  // index of the first basis function
  inline int  operator() (double x, TinyVector<double,4>& bfuncs) const;
  // Same as above, but also computes first derivatives
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv) const;
  // Same as above, but also computes second derivatives
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv,
			  TinyVector<double,4> &deriv2) const;
  inline int  operator() (double x, TinyVector<float,4>& bfunc,
			  TinyVector<float,4> &deriv,
			  TinyVector<float,4> &deriv2) const;
  // These versions take the grid point index.  These are needed for
  // the boundary conditions on the interpolating equations.  The
  // complex version return the basis functions in both the real and
  // imaginary parts.
  inline void operator() (int i, TinyVector<float,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs) const;

  inline void operator() (int i, TinyVector<float,4>& bfuncs,
			  TinyVector<float,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs) const;

  inline void operator() (int i, TinyVector<float,4>& bfuncs,
			  TinyVector<float,4> &dbfuncs,
			  TinyVector<float,4> &d2bfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs,
			  TinyVector<double,4> &d2bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs,
			  TinyVector<complex<double>,4> &d2bfuncs) const;
  inline NUBsplineBasis();
};

#ifdef __SSE3__
#include <xmmintrin.h>
inline int
NUBsplineBasis<LinearGrid>::Find(double x) const
{
  double delta = (x - GridStart);
  double di = Linv * delta;
  // Compute floor of di with SSE instructions
  int floor_di;
  __m128d _di, _floor, _zero, _mask, _m1;
  _di = _mm_set_sd (di);
  _m1 = _mm_set_sd (-1.0);
  _zero = _mm_setzero_pd();
  _mask = _mm_cmplt_sd (_di, _zero);
  _m1   = _mm_and_pd(_m1, _mask);
  _di   = _mm_add_sd (_di, _m1);
  floor_di =_mm_cvttsd_si32 (_di);
  //  cerr << "di = " << di << "  floor(di) = " << floor_di << endl;
  delta -= Periodic * (double)floor_di*L;
  double fi = delta * GridDeltaInv;
  __m128d _fi = _mm_set_sd(fi);
  int i = _mm_cvttsd_si32 (_fi);
  double t = fi - (double)i;
  tp[3] = 1.0;
  tp[2] = t;
  tp[1] = t*t;
  tp[0] = t*t*t;
  _tp = _mm_set_ps (tp[0], tp[1], tp[2], tp[3]);
  return i;
}

#else

inline int
NUBsplineBasis<LinearGrid>::Find(double x) const
{
  double delta = x - GridStart;
  delta -= Periodic*floor(delta*Linv)*L;
  double fi = delta*GridDeltaInv;
  double ipart;
  double t = modf (fi, &ipart);
  int i = (int) ipart;
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  return i;
}

#endif

inline int
NUBsplineBasis<LinearGrid>::operator()(double x, TinyVector<double,4> &bfuncs) const
{
  int i0 = Find(x);
  bfuncs[0]  = A(0,0)*tp[0] + A(1,0)*tp[1] + A(2,0)*tp[2] + A(3,0)*tp[3];
  bfuncs[1]  = A(0,1)*tp[0] + A(1,1)*tp[1] + A(2,1)*tp[2] + A(3,1)*tp[3];
  bfuncs[2]  = A(0,2)*tp[0] + A(1,2)*tp[1] + A(2,2)*tp[2] + A(3,2)*tp[3];
  bfuncs[3]  = A(0,3)*tp[0] + A(1,3)*tp[1] + A(2,3)*tp[2] + A(3,3)*tp[3];

  return i0;
}



inline int
NUBsplineBasis<LinearGrid>::operator()(double x, TinyVector<double,4> &bfuncs,
				       TinyVector<double,4> &dbfuncs) const
{
  int i0 = Find(x);
  bfuncs[0]  = A(0,0)*tp[0] + A(1,0)*tp[1] + A(2,0)*tp[2] + A(3,0)*tp[3];
  bfuncs[1]  = A(0,1)*tp[0] + A(1,1)*tp[1] + A(2,1)*tp[2] + A(3,1)*tp[3];
  bfuncs[2]  = A(0,2)*tp[0] + A(1,2)*tp[1] + A(2,2)*tp[2] + A(3,2)*tp[3];
  bfuncs[3]  = A(0,3)*tp[0] + A(1,3)*tp[1] + A(2,3)*tp[2] + A(3,3)*tp[3];

  dbfuncs[0]  = GridDeltaInv*(dA(1,0)*tp[1] + dA(2,0)*tp[2] + dA(3,0)*tp[3]);
  dbfuncs[1]  = GridDeltaInv*(dA(1,1)*tp[1] + dA(2,1)*tp[2] + dA(3,1)*tp[3]);
  dbfuncs[2]  = GridDeltaInv*(dA(1,2)*tp[1] + dA(2,2)*tp[2] + dA(3,2)*tp[3]);
  dbfuncs[3]  = GridDeltaInv*(dA(1,3)*tp[1] + dA(2,3)*tp[2] + dA(3,3)*tp[3]);

  return i0;
}

inline int
NUBsplineBasis<LinearGrid>::operator()(double x, TinyVector<double,4> &bfuncs,
				       TinyVector<double,4> &dbfuncs,
				       TinyVector<double,4> &d2bfuncs) const
{
  int i0 = Find(x);
  bfuncs[0]  = A(0,0)*tp[0] + A(1,0)*tp[1] + A(2,0)*tp[2] + A(3,0)*tp[3];
  bfuncs[1]  = A(0,1)*tp[0] + A(1,1)*tp[1] + A(2,1)*tp[2] + A(3,1)*tp[3];
  bfuncs[2]  = A(0,2)*tp[0] + A(1,2)*tp[1] + A(2,2)*tp[2] + A(3,2)*tp[3];
  bfuncs[3]  = A(0,3)*tp[0] + A(1,3)*tp[1] + A(2,3)*tp[2] + A(3,3)*tp[3];

  dbfuncs[0]  = GridDeltaInv*(dA(1,0)*tp[1] + dA(2,0)*tp[2] + dA(3,0)*tp[3]);
  dbfuncs[1]  = GridDeltaInv*(dA(1,1)*tp[1] + dA(2,1)*tp[2] + dA(3,1)*tp[3]);
  dbfuncs[2]  = GridDeltaInv*(dA(1,2)*tp[1] + dA(2,2)*tp[2] + dA(3,2)*tp[3]);
  dbfuncs[3]  = GridDeltaInv*(dA(1,3)*tp[1] + dA(2,3)*tp[2] + dA(3,3)*tp[3]);

  d2bfuncs[0]  = GridDeltaInv*GridDeltaInv*(d2A(2,0)*tp[2] + d2A(3,0)*tp[3]);
  d2bfuncs[1]  = GridDeltaInv*GridDeltaInv*(d2A(2,1)*tp[2] + d2A(3,1)*tp[3]);
  d2bfuncs[2]  = GridDeltaInv*GridDeltaInv*(d2A(2,2)*tp[2] + d2A(3,2)*tp[3]);
  d2bfuncs[3]  = GridDeltaInv*GridDeltaInv*(d2A(2,3)*tp[2] + d2A(3,3)*tp[3]);
  return i0;
}

#ifdef __SSE3__
inline int
NUBsplineBasis<LinearGrid>::operator()(double x, TinyVector<float,4> &bfuncs,
				       TinyVector<float,4> &dbfuncs,
				       TinyVector<float,4> &d2bfuncs) const
{
  int i0 = Find(x);
  __m128 tmp0, tmp1, tmp2, tmp3, GDE;
  tmp0 = _mm_mul_ps (_A[0], _tp);
  tmp1 = _mm_mul_ps (_A[1], _tp);
  tmp2 = _mm_mul_ps (_A[2], _tp);
  tmp3 = _mm_mul_ps (_A[3], _tp);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp1 = _mm_hadd_ps (tmp2, tmp3);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  _mm_store_ps (&(bfuncs[0]), tmp0);

  tmp0 = _mm_mul_ps (_dA[0], _tp);
  tmp1 = _mm_mul_ps (_dA[1], _tp);
  tmp2 = _mm_mul_ps (_dA[2], _tp);
  tmp3 = _mm_mul_ps (_dA[3], _tp);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp1 = _mm_hadd_ps (tmp2, tmp3);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp0 = _mm_mul_ps  (tmp0, _GDE);
  _mm_store_ps (&(dbfuncs[0]), tmp0);

  tmp0 = _mm_mul_ps (_d2A[0], _tp);
  tmp1 = _mm_mul_ps (_d2A[1], _tp);
  tmp2 = _mm_mul_ps (_d2A[2], _tp);
  tmp3 = _mm_mul_ps (_d2A[3], _tp);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp1 = _mm_hadd_ps (tmp2, tmp3);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp0 = _mm_mul_ps  (tmp0, _GDE);
  tmp0 = _mm_mul_ps  (tmp0, _GDE);
  _mm_store_ps (&(d2bfuncs[0]), tmp0);

  return i0;
}

#else
inline int
NUBsplineBasis<LinearGrid>::operator()(double x, TinyVector<float,4> &bfuncs,
				       TinyVector<float,4> &dbfuncs,
				       TinyVector<float,4> &d2bfuncs) const
{
  


}


#endif


inline void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<double,4> &bfuncs) const
{
  bfuncs = TinyVector<double,4> (Sixth, TwoThirds, Sixth, 0.0);
}



inline void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<double,4> &bfuncs,
				       TinyVector<double,4> &dbfuncs) const
{
  bfuncs  = TinyVector<double,4> (Sixth, TwoThirds, Sixth, 0.0);
  dbfuncs = TinyVector<double,4> (-0.5, 0.0, 0.5, 0.0);
}

inline void
NUBsplineBasis<LinearGrid>::operator()(int i,
				       TinyVector<complex<double>,4> &bfuncs) const
{
  TinyVector<double,4> funcs;
  (*this)(i, funcs);
  bfuncs[0] = complex<double>(funcs[0], funcs[0]);
  bfuncs[1] = complex<double>(funcs[1], funcs[1]);
  bfuncs[2] = complex<double>(funcs[2], funcs[2]);
  bfuncs[3] = complex<double>(funcs[3], funcs[3]);
}

inline void
NUBsplineBasis<LinearGrid>::operator()(int i,
				       TinyVector<float,4> &bfuncs) const
{
  TinyVector<double,4> funcs;
  (*this)(i, funcs);
  bfuncs[0] = float(funcs[0]);
  bfuncs[1] = float(funcs[1]);
  bfuncs[2] = float(funcs[2]);
  bfuncs[3] = float(funcs[3]);
}



void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				     TinyVector<complex<double>,4> &dbfuncs) const
{
  TinyVector<double,4> funcs, dfuncs;
  (*this)(i, funcs, dfuncs);
  bfuncs[0]  = complex<double>(funcs[0],  funcs[0]);
  bfuncs[1]  = complex<double>(funcs[1],  funcs[1]);
  bfuncs[2]  = complex<double>(funcs[2],  funcs[2]);
  bfuncs[3]  = complex<double>(funcs[3],  funcs[3]);
  dbfuncs[0] = complex<double>(dfuncs[0], dfuncs[0]);
  dbfuncs[1] = complex<double>(dfuncs[1], dfuncs[1]);
  dbfuncs[2] = complex<double>(dfuncs[2], dfuncs[2]);
  dbfuncs[3] = complex<double>(dfuncs[3], dfuncs[3]);
}

#ifdef __SSE3__
void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<float,4> &bfuncs,
				     TinyVector<float,4> &dbfuncs) const
{
  __m128 tmp0, tmp1, tmp2, tmp3;
  tmp0 = _mm_mul_ps (_A[0], _tp);
  tmp1 = _mm_mul_ps (_A[1], _tp);
  tmp2 = _mm_mul_ps (_A[2], _tp);
  tmp3 = _mm_mul_ps (_A[3], _tp);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  tmp1 = _mm_hadd_ps (tmp2, tmp3);
  tmp0 = _mm_hadd_ps (tmp0, tmp1);
  _mm_store_ps (&(bfuncs[0]), tmp0);
}

#else
void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<float,4> &bfuncs,
				     TinyVector<float,4> &dbfuncs) const
{
  TinyVector<double,4> funcs, dfuncs;
  (*this)(i, funcs, dfuncs);
  bfuncs[0]  = float(funcs[0]);
  bfuncs[1]  = float(funcs[1]);
  bfuncs[2]  = float(funcs[2]);
  bfuncs[3]  = float(funcs[3]);
  dbfuncs[0] = float(dfuncs[0]);
  dbfuncs[1] = float(dfuncs[1]);
  dbfuncs[2] = float(dfuncs[2]);
  dbfuncs[3] = float(dfuncs[3]);
}
#endif


inline void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<double,4> &bfuncs,
				       TinyVector<double,4> &dbfuncs,
				       TinyVector<double,4> &d2bfuncs) const
{
  bfuncs   = TinyVector<double,4> (Sixth, TwoThirds, Sixth, 0.0);
  dbfuncs  = TinyVector<double,4> (-0.5, 0.0, 0.5, 0.0);
  d2bfuncs = TinyVector<double,4> (1.0, -2.0, 1.0, 0.0);
}


inline void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				       TinyVector<complex<double>,4> &dbfuncs,
				       TinyVector<complex<double>,4> &d2bfuncs) const
{
  TinyVector<double,4> funcs, dfuncs, d2funcs;
  (*this)(i, funcs, dfuncs, d2funcs);
  bfuncs[0]   = complex<double>(  funcs[0],   funcs[0]);
  bfuncs[1]   = complex<double>(  funcs[1],   funcs[1]);
  bfuncs[2]   = complex<double>(  funcs[2],   funcs[2]);
  bfuncs[3]   = complex<double>(  funcs[3],   funcs[3]);
  dbfuncs[0]  = complex<double>( dfuncs[0],  dfuncs[0]);
  dbfuncs[1]  = complex<double>( dfuncs[1],  dfuncs[1]);
  dbfuncs[2]  = complex<double>( dfuncs[2],  dfuncs[2]);
  dbfuncs[3]  = complex<double>( dfuncs[3],  dfuncs[3]);
  d2bfuncs[0] = complex<double>(d2funcs[0], d2funcs[0]);
  d2bfuncs[1] = complex<double>(d2funcs[1], d2funcs[1]);
  d2bfuncs[2] = complex<double>(d2funcs[2], d2funcs[2]);
  d2bfuncs[3] = complex<double>(d2funcs[3], d2funcs[3]);
}

inline void
NUBsplineBasis<LinearGrid>::operator()(int i, TinyVector<float,4> &bfuncs,
				       TinyVector<float,4> &dbfuncs,
				       TinyVector<float,4> &d2bfuncs) const
{
  TinyVector<double,4> funcs, dfuncs, d2funcs;
  (*this)(i, funcs, dfuncs, d2funcs);
  bfuncs[0]   = float(  funcs[0]);
  bfuncs[1]   = float(  funcs[1]);
  bfuncs[2]   = float(  funcs[2]);
  bfuncs[3]   = float(  funcs[3]);
  dbfuncs[0]  = float( dfuncs[0]);
  dbfuncs[1]  = float( dfuncs[1]);
  dbfuncs[2]  = float( dfuncs[2]);
  dbfuncs[3]  = float( dfuncs[3]);
  d2bfuncs[0] = float(d2funcs[0]);
  d2bfuncs[1] = float(d2funcs[1]);
  d2bfuncs[2] = float(d2funcs[2]);
  d2bfuncs[3] = float(d2funcs[3]);
}


inline void
NUBsplineBasis<LinearGrid>::Init(LinearGrid *grid, bool periodic)
{
  Periodic = periodic ? 1.0 : 0.0;
  GridStart = grid->Start;
  GridEnd   = grid->End;
  int N = grid->NumPoints;
  L = GridEnd - GridStart;
  Linv = 1.0/L;
  GridDelta = L/(double)(N-1);
  GridDeltaInv = 1.0/GridDelta;
#ifdef __SSE3__
  _GDE = _mm_set_ps (GridDeltaInv, GridDeltaInv, GridDeltaInv, GridDeltaInv);
#endif
}


NUBsplineBasis<LinearGrid>::NUBsplineBasis()
{
  Sixth     = 1.0/6.0;
  TwoThirds = 2.0/3.0;
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;

  dA(0,0)= 0.0; dA(0,1)= 0.0; dA(0,2)= 0.0; dA(0,3)= 0.0;
  dA(1,0)=-0.5; dA(1,1)= 1.5; dA(1,2)=-1.5; dA(1,3)= 0.5;
  dA(2,0)= 1.0; dA(2,1)=-2.0; dA(2,2)= 1.0; dA(2,3)= 0.0;
  dA(3,0)=-0.5; dA(3,1)= 0.0; dA(3,2)= 0.5; dA(3,3)= 0.0;

  d2A(0,0)= 0.0; d2A(0,1)= 0.0; d2A(0,2)= 0.0; d2A(0,3)= 0.0;
  d2A(1,0)= 0.0; d2A(1,1)= 0.0; d2A(1,2)= 0.0; d2A(1,3)= 0.0;
  d2A(2,0)=-1.0; d2A(2,1)= 3.0; d2A(2,2)=-3.0; d2A(2,3)= 1.0;
  d2A(3,0)= 1.0; d2A(3,1)=-2.0; d2A(3,2)= 1.0; d2A(3,3)= 0.0;

  d3A(0,0)= 0.0; d3A(0,1)= 0.0; d3A(0,2)= 0.0; d3A(0,3)= 0.0;
  d3A(1,0)= 0.0; d3A(1,1)= 0.0; d3A(1,2)= 0.0; d3A(1,3)= 0.0;
  d3A(2,0)= 0.0; d3A(2,1)= 0.0; d3A(1,2)= 2.0; d3A(2,3)= 0.0;
  d3A(3,0)=-1.0; d3A(3,1)= 3.0; d3A(3,2)=-3.0; d3A(3,3)= 1.0;
#ifdef __SSE3__
  _A[0] = _mm_set_ps (-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0);
  _A[1] = _mm_set_ps ( 3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0);
  _A[2] = _mm_set_ps (-3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0);
  _A[3] = _mm_set_ps ( 1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0);

  _dA[0] = _mm_set_ps ( 0.0, -0.5,  1.0, -0.5);
  _dA[1] = _mm_set_ps ( 0.0,  1.5, -2.0,  0.0);
  _dA[2] = _mm_set_ps ( 0.0, -1.5,  1.0,  0.5);
  _dA[3] = _mm_set_ps ( 0.0,  0.5,  0.0,  0.0);

  _d2A[0] = _mm_set_ps ( 0.0,  0.0, -1.0,  1.0);
  _d2A[1] = _mm_set_ps ( 0.0,  0.0,  3.0, -2.0);
  _d2A[2] = _mm_set_ps ( 0.0,  0.0, -3.0,  1.0);
  _d2A[3] = _mm_set_ps ( 0.0,  0.0,  1.0,  0.0);

#endif
}



template<typename GridType, typename T> inline void
SolveDerivInterp1D (NUBsplineBasis<GridType> &basis,
		    Array<T,1> data, Array<T,1> p,
		    TinyVector<T,4> abcdInitial,
		    TinyVector<T,4> abcdFinal)
{
  assert (p.size() == (data.size()+2));

  // Banded matrix storage.  The first three elements in the
  // tinyvector store the tridiagonal coefficients.  The last element
  // stores the RHS data.
  Array<TinyVector<T,4>,1> bands(p.size());
  int M = data.size();

  // Fill up bands
  bands(0)          = abcdInitial;
  bands(p.size()-1) = abcdFinal;
  for (int i=0; i<M; i++) {
    basis (i, bands(i+1));
    bands(i+1)[3] = data(i);
  }

  // Now solve:
  // First and last rows are different
  bands(0)[1] /= bands(0)[0];
  bands(0)[2] /= bands(0)[0];
  bands(0)[3] /= bands(0)[0];
  bands(0)[0] = 1.0;
  bands(1)[1] -= bands(1)[0]*bands(0)[1];
  bands(1)[2] -= bands(1)[0]*bands(0)[2];
  bands(1)[3] -= bands(1)[0]*bands(0)[3];
  bands(0)[0] = 0.0;
  bands(1)[2] /= bands(1)[1];
  bands(1)[3] /= bands(1)[1];
  bands(1)[1] = 1.0;
  
  // Now do rows 2 through M+1
  for (int row=2; row < (p.size()-1); row++) {
    bands(row)[1] -= bands(row)[0]*bands(row-1)[2];
    bands(row)[3] -= bands(row)[0]*bands(row-1)[3];
    bands(row)[2] /= bands(row)[1];
    bands(row)[3] /= bands(row)[1];
    bands(row)[0] = 0.0;
    bands(row)[1] = 1.0;
  }

  // Do last row
  bands(M+1)[1] -= bands(M+1)[0]*bands(M-1)[2];
  bands(M+1)[3] -= bands(M+1)[0]*bands(M-1)[3];
  bands(M+1)[2] -= bands(M+1)[1]*bands(M)[2];
  bands(M+1)[3] -= bands(M+1)[1]*bands(M)[3];
  bands(M+1)[3] /= bands(M+1)[2];
  bands(M+1)[2] = 1.0;

  p(M+1) = bands(M+1)[3];
  // Now back substitute up
  for (int row=M; row>0; row--)
    p(row) = bands(row)[3] - bands(row)[2]*p(row+1);
  
  // Finish with first row
  p(0) = bands(0)[3] - bands(0)[1]*p(1) - bands(0)[2]*p(2);
}

template<typename T, int N> inline TinyVector<T,4> 
real(TinyVector<complex<T>,N> v)
{
  TinyVector<T,N> rv;
  for (int i=0; i<N; i++)
    rv[i] = real(v[i]);
  return rv;
}

template<typename T, int N> inline TinyVector<T,4> 
imag(TinyVector<complex<T>,N> v)
{
  TinyVector<T,N> iv;
  for (int i=0; i<N; i++)
    iv[i] = imag(v[i]);
  return iv;
}


template<typename GridType, typename T> inline void
SolveDerivInterp1D (NUBsplineBasis<GridType> &basis,
		    Array<complex<T>,1> data, Array<complex<T>,1> p,
		    TinyVector<complex<T>,4> abcdInitial,
		    TinyVector<complex<T>,4> abcdFinal)
{
  Array<T,1> rdata(data.size()), rp(p.size()),
    idata(data.size()), ip(p.size());
  TinyVector<T,4> rStartBC, iStartBC, rEndBC, iEndBC;
  rStartBC = real(abcdInitial);
  iStartBC = imag(abcdInitial);
  rEndBC   = real(abcdFinal);
  iEndBC   = imag(abcdFinal);
  
  for (int i=0; i<data.size(); i++){
    rdata(i) = real(data(i));
    idata(i) = imag(data(i));
  }
  SolveDerivInterp1D (basis, rdata, rp, rStartBC, rEndBC);
  SolveDerivInterp1D (basis, idata, ip, iStartBC, iEndBC);
  for (int i=0; i<p.size(); i++)
    p(i) = complex<T> (rp(i), ip(i));
}


template<typename GridType, typename T> inline void
SolvePeriodicInterp1D (NUBsplineBasis<GridType> &basis,
		       Array<T,1> data, Array<T,1> p)
{
  assert (p.size() == (data.size()+3));

  // Banded matrix storage.  The first three elements in the
  // tinyvector store the tridiagonal coefficients.  The last element
  // stores the RHS data.
  Array<TinyVector<double,4>,1> bands(data.size());
  Array<double,1> lastCol(data.size());
  int M = data.size();

  // Fill up bands
  for (int i=0; i<M; i++) {
    basis (i, bands(i));
    bands(i)[3] = data(i);
  }
    
  // Now solve:
  // First and last rows are different
  bands(0)[2] /= bands(0)[1];
  bands(0)[0] /= bands(0)[1];
  bands(0)[3] /= bands(0)[1];
  bands(0)[1]  = 1.0;
  bands(M-1)[1] -= bands(M-1)[2]*bands(0)[0];
  bands(M-1)[3] -= bands(M-1)[2]*bands(0)[3];
  bands(M-1)[2]  = -bands(M-1)[2]*bands(0)[2];
  lastCol(0) = bands(0)[0];
  
  for (int row=1; row < (M-1); row++) {
    bands(row)[1] -= bands(row)[0] * bands(row-1)[2];
    bands(row)[3] -= bands(row)[0] * bands(row-1)[3];
    lastCol(row)   = -bands(row)[0] * lastCol(row-1);
    bands(row)[0] = 0.0;
    bands(row)[2] /= bands(row)[1];
    bands(row)[3] /= bands(row)[1];
    lastCol(row)  /= bands(row)[1];
    bands(row)[1]  = 1.0;
    if (row < (M-2)) {
      bands(M-1)[3] -= bands(M-1)[2]*bands(row)[3];
      bands(M-1)[1] -= bands(M-1)[2]*lastCol(row);
      bands(M-1)[2] = -bands(M-1)[2]*bands(row)[2];
    }
  }

  // Now do last row
  // The [2] element and [0] element are now on top of each other 
  bands(M-1)[0] += bands(M-1)[2];
  bands(M-1)[1] -= bands(M-1)[0] * (bands(M-2)[2]+lastCol(M-2));
  bands(M-1)[3] -= bands(M-1)[0] *  bands(M-2)[3];
  bands(M-1)[3] /= bands(M-1)[1];
  p(M) = bands(M-1)[3];
  for (int row=M-2; row>=0; row--) 
    p(row+1) = bands(row)[3] - bands(row)[2]*p(row+2) - lastCol(row)*p(M);
  
  p(0) = p(M);
  p(M+1) = p(1);
  p(M+2) = p(2);
}


template<typename GridType, typename T> inline void
SolvePeriodicInterp1D (NUBsplineBasis<GridType> &basis,
		       Array<complex<T>,1> data, Array<complex<T>,1> p)
{
  Array<T,1> rdata(data.size()), rp(p.size()),
    idata(data.size()), ip(p.size());
  
  for (int i=0; i<data.size(); i++){
    rdata(i) = real(data(i));
    idata(i) = imag(data(i));
  }
  SolvePeriodicInterp1D (basis, rdata, rp);
  SolvePeriodicInterp1D (basis, idata, ip);
  for (int i=0; i<p.size(); i++)
    p(i) = complex<T> (rp(i), ip(i));
}



#endif
