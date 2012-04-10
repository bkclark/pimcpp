#include "TricubicBspline.h"

// Solve for fixed first derivative.  
// data(0)  = dy_0/dx dx.  
// data(1)   = y_0 and
// data(M-1) = y_{M-1}.  
// data(M)   = dy_{M-1}/dx dx.  
template<typename T> inline void
SolveFirstDerivInterp1D (Array<T,1> data, Array<T,1> p)
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
  mu(1) = -1.0;
  d(0) /= (-ratio2);
  for (int row=1; row <=M; row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    mu(row) *= diagInv;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
  }
  d(M+1) -= -ratio2 * data(M-1);
  mu(M+1) -= -ratio2*mu(M-1);
  double diag = ratio2 - mu(M+1)*mu(M);
  p(M+1) = d(M+1)/diag;
 
  // Now go back upward, back substituting
  for (int row=M; row>=1; row--) 
    p(row+1) = d(row) - mu(row)*p(row+1) - gamma(row)*p(M-1);

  // And do 0th row
  p(0) = d(0) + p(2)*d(2);
}



template<typename T> void
TricubicBspline<T>::SolvePeriodicInterp (const Array<T,3> &data)
{
  // Do X direction
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) 
      SolvePeriodicInterp1D(data(Range(0,Nx-1), iy, iz), P(Range(1,Nx),iy+1, iz+1));
  
  // Do Y direction
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) 
      SolvePeriodicInterp1D(P(ix+1,Range(1,Ny), iz+1), P(ix+1, Range(1,Ny), iz+1));
  
  // Do z direction
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) 
      SolvePeriodicInterp1D(P(ix+1,iy+1,Range(1,Nz)), P(ix+1, iy+1, Range(1,Nz)));
}

// template<typename T> void
// TricubicBspline<T>::SolveClampedInterp (Array<T,3> &data)
// {
//   Array<T,1> dx(data.extent(0)+2);
//   Array<T,1> dy(data.extent(0)+2);
//   // Do X direction
//   for (int iy=0; iy<Ny; iy++)
//     for (int iz=0; iz<Nz; iz++) 
//       SolvePeriodicInterp1D(data(Range(0,Nx-1), iy, iz), P(Range(1,Nx),iy+1, iz+1));
  
//   // Do Y direction
//   for (int ix=0; ix<Nx; ix++)
//     for (int iz=0; iz<Nz; iz++) 
//       SolvePeriodicInterp1D(P(ix+1,Range(1,Ny), iz+1), P(ix+1, Range(1,Ny), iz+1));
  
//   // Do z direction
//   for (int ix=0; ix<Nx; ix++)
//     for (int iy=0; iy<Ny; iy++) 
//       SolvePeriodicInterp1D(P(ix+1,iy+1,Range(1,Nz)), P(ix+1, iy+1, Range(1,Nz)));
// }

template<typename T> void
TricubicBspline<T>::SolveInterp (const Array<T,3> &data,
				 BCType xbc, BCType ybc, BCType zbc)
{
  int Mx, My, Mz;
  Mx = (xbc == PERIODIC) ? Nx+3 : Nx+2;
  My = (ybc == PERIODIC) ? Ny+3 : Ny+2;
  Mz = (zbc == PERIODIC) ? Nz+3 : Nz+2;
  
  dx = Lx/(double)(Mx-3); dxInv = 1.0/dx;
  dy = Ly/(double)(My-3); dyInv = 1.0/dy;
  dz = Lz/(double)(Mz-3); dzInv = 1.0/dz;

  P.resize(Mx, My, Mz);

  TinyVector<double,4> BC;
  ////////////////////
  // Do X direction //
  ////////////////////
  if (xbc == PERIODIC) {
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	SolvePeriodicInterp1D(data(Range::all(), iy, iz), P(Range::all(),iy+1, iz+1));
  }
  else {
    if (xbc == FLAT) 
      BC = -3.0, 0.0, 3.0, 0.0;
    else if (xbc == NATURAL)
      BC = 6.0, -12.0, 6.0, 0.0;
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	SolveDerivInterp1D(data(Range::all(), iy, iz),  P(Range::all(), iy+1, iz+1), BC, BC);
  }
  ////////////////////
  // Do Y direction //
  ////////////////////
  if (ybc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++) 
	SolvePeriodicInterp1D(P(ix, Range(1,Ny), iz), P(ix, Range::all(), iz));
  }
  else {
    if (ybc == FLAT) 
      BC = -3.0, 0.0, 3.0, 0.0;
    else if (ybc == NATURAL)
      BC = 6.0, -12.0, 6.0, 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++)
	SolveDerivInterp1D(P(ix,Range(1,Ny), iz),  P(ix, Range::all(), iz), BC, BC);
  }
  ////////////////////
  // Do Z direction //
  ////////////////////
  if (zbc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++) 
	SolvePeriodicInterp1D(P(ix, iy, Range(1,Nz)), P(ix, iy, Range::all()));
  }
  else {
    if (zbc == FLAT) 
      BC = -3.0, 0.0, 3.0, 0.0;
    else if (zbc == NATURAL)
      BC = 6.0, -12.0, 6.0, 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++)
	SolveDerivInterp1D(P(ix, iy, Range(1,Nz)),  P(ix, iy, Range::all()), BC, BC);
  }
}


template<typename T> void
TricubicBspline<T>::MakePeriodic()
{
  // Now, make periodic
  for (int ix=0; ix<(Nx+3); ix++)
    for (int iy=0; iy<(Ny+3); iy++)
      for (int iz=0; iz<(Nz+3); iz++)
	P(ix, iy, iz) = P((ix+Nx-1)%Nx+1, (iy+Ny-1)%Ny+1, (iz+Nz-1)%Nz+1);
}


template<typename T> void
TricubicBspline<T>::Init (double xi, double xf, double yi, double yf, double zi, double zf,
			  const Array<T,3> &data, bool interp, 
			  BCType xbc, BCType ybc, BCType zbc)
{
  if ((xbc==FIXED_FIRST) || (xbc==FIXED_SECOND) ||
      (ybc==FIXED_FIRST) || (ybc==FIXED_SECOND) ||
      (zbc==FIXED_FIRST) || (zbc==FIXED_SECOND)) {
    cerr << "Cannot fix nonzero first or second derivative wint TricubicBSpline without specifying.\n";
    abort();
  }
  Nx = data.extent(0);
  Ny = data.extent(1);
  Nz = data.extent(2);
  xStart=xi;  xEnd=xf;  yStart=yi;  yEnd=yf;  zStart=zi; zEnd=zf;
  Lx    = xf-xi;  Ly    = yf-yi;  Lz    = zf-zi;
  LxInv = 1.0/Lx; LyInv = 1.0/Ly; LzInv = 1.0/Lz;
//   dx = Lx/(double)Nx; dxInv = 1.0/dx;
//   dy = Ly/(double)Ny; dyInv = 1.0/dy;
//   dz = Lz/(double)Nz; dzInv = 1.0/dz;
//   P.resize(Nx+3, Ny+3, Nz+3);

  Interpolating = interp;

  if (interp) 
    SolveInterp (data, xbc, ybc, zbc);
  else if ((xbc==PERIODIC) && (ybc==PERIODIC) && (zbc==PERIODIC)) {
    P.resize(Nx+3, Ny+3, Nz+3);
    dx = Lx/(double)Nx; dxInv = 1.0/dx;
    dy = Ly/(double)Ny; dyInv = 1.0/dy;
    dz = Lz/(double)Nz; dzInv = 1.0/dz;
    P(Range(1,Nx),Range(1,Ny),Range(1,Nz)) = data;
    MakePeriodic();
  }
  else {
    cerr << "I don't know how to construct a nonperiodic, noninterpolating Tricubic B-spline.\n";
    abort();
  }

//   if (xbc==PERIODIC) {
//     if (interp)
//       SolvePeriodicInterp(data);
//     else
//       P(Range(1,Nx),Range(1,Ny),Range(1,Nz)) = data;
//     MakePeriodic();
//   }
//   else {
//     cerr << "Nonperiodic Tricubic B-splines not yet supported.\n";
//     abort();
//   }
}

template<typename T>
TricubicBspline<T>::TricubicBspline()
{
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
}


template class TricubicBspline<double>;
template class TricubicBspline<float>;
template class TricubicBspline<complex<double> >;

