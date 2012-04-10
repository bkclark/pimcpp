#ifndef CUBIC_NUB_SPLINE_H
#define CUBIC_NUB_SPLINE_H

#include "NUBsplineBasis.h"

template<typename GridType>
class CubicNUBspline
{
private:
  NUBsplineBasis<GridType> Basis;
  bool Interpolating, Periodic;
  GridType *SplineGrid;

  // The control points
  Array<double,1> P;

public:
  inline double GetControlPoint (int i);
  void Init(GridType &grid, Array<double,1> &data, 
	    BoundaryCondition<double> startBC=BoundaryCondition<double>(PERIODIC),
	    BoundaryCondition<double>   endBC=BoundaryCondition<double>(PERIODIC));
  inline double operator()(double x) const;
  inline double Deriv     (double x) const;
  inline double Deriv2    (double x) const;
  inline double Deriv3    (double x) const;

  CubicNUBspline() : Periodic(false), Interpolating(false)
  { }
};


template<typename GridType> void
CubicNUBspline<GridType>::Init(GridType &grid, Array<double,1> &data,
			       BoundaryCondition<double> startBC,
			       BoundaryCondition<double>   endBC)
{
  TinyVector<double,4> rBC, lBC, dummy1, dummy2;
  int N = grid.NumPoints;
  SplineGrid = &grid;

  if (startBC.GetType() == PERIODIC) {
    Basis.Init(&grid, true);
  
    assert (endBC.GetType() == PERIODIC);
    P.resize(data.size()+3);
    assert (grid.NumPoints == data.size()+1);
    
    SolvePeriodicInterp1D (Basis, data, P);
  }
  else {
    Basis.Init(&grid);
    if (startBC.GetType() == FIXED_FIRST) 
      Basis(0, dummy1, lBC);
    else if (startBC.GetType() == FIXED_SECOND) 
      Basis(0, dummy1, dummy2, lBC);
    else {
      cerr << "Unknown right BC type.\n";
      abort();
    }
    lBC[3] = startBC.GetVal();
    
    if (endBC.GetType() == FIXED_FIRST) 
      Basis(N-1, dummy1, rBC);
    else if (endBC.GetType() == FIXED_SECOND) 
      Basis(N-1, dummy1, dummy2, rBC);
    else {
      cerr << "Unknown right BC type.\n";
      abort();
    }
    rBC[3] = endBC.GetVal();
    
    P.resize(data.size()+2);
    assert (grid.NumPoints == data.size());
    
    SolveDerivInterp1D (Basis, data, P, lBC, rBC);
  }
}

template<typename GridType>
inline double
CubicNUBspline<GridType>::operator()(double x) const
{
  TinyVector<double,4> bfuncs;
  x = min (x, SplineGrid->End - 1.0e-10);
  int i = Basis (x, bfuncs);
  return P(i)*bfuncs[0] + P(i+1)*bfuncs[1] + P(i+2)*bfuncs[2] + P(i+3)*bfuncs[3];
}

template<typename GridType>
inline double
CubicNUBspline<GridType>::Deriv(double x) const
{
  TinyVector<double,4> bfuncs, dbfuncs;
  int i = Basis (x, bfuncs, dbfuncs);
  return (P(i)*dbfuncs[0] + P(i+1)*dbfuncs[1] + 
	  P(i+2)*dbfuncs[2] + P(i+3)*dbfuncs[3]);
}

template<typename GridType>
inline double
CubicNUBspline<GridType>::Deriv2(double x) const
{
  TinyVector<double,4> bfuncs, dbfuncs, d2bfuncs;
  int i = Basis (x, bfuncs, dbfuncs, d2bfuncs);
  return (P(i)*d2bfuncs[0] + P(i+1)*d2bfuncs[1] + 
	  P(i+2)*d2bfuncs[2] + P(i+3)*d2bfuncs[3]);
}



#endif
