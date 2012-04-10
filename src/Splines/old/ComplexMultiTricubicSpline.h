#ifndef COMPLEX_MULTI_TRICUBIC_SPLINE_H
#define COMPLEX_MULTI_TRICUBIC_SPLINE_H

#include "Grid.h"
#include <cmath>
//#include <blitz/array.h>
//using namespace blitz;

/// Each point of F contains:
/// 0)  real [F(x,y,z)]
/// 1)  real [dF/dx]
/// 2)  real [dF/dy]
/// 3)  real [dF/dz]
/// 4)  real [d2F/dxdy]
/// 5)  real [d2F/dxdz]
/// 6)  real [d2F/dydz]
/// 7)  real [d3F/dxdydz]
/// 8)  imag [F(x,y,z)]
/// 9)  imag [dF/dx]
/// 10) imag [dF/dy]
/// 11) imag [dF/dz]
/// 12) imag [d2F/dxdy]
/// 13) imag [d2F/dxdz]
/// 14) imag [d2F/dydz]
/// 15) imag [d3F/dxdydz]


class ComplexMultiTricubicSpline
{
  inline double p1(double t)
  { return ((t-1.0)*(t-1.0)*(1.0+2.0*t)); }
  inline double p2(double t)
  { return (t*t*(3.0-2.0*t)); }
  inline double q1(double t)
  { return (t*(t-1.0)*(t-1.0)); }
  inline double q2(double t)
  { return (t*t*(t-1.0)); }
  inline double dp1(double t)
  { return (6.0*t*(t-1.0)); }
  inline double dq1(double t)
  { return ((t-1.0)*(3.0*t-1.0)); }
  inline double dp2(double t)
  { return (-dp1(t)); }
  inline double dq2 (double t)
  { return (3.0*t*t - 2.0*t); }
  inline double d2p1(double t)
  { return (12.0*t-6.0); }
  inline double d2q1 (double t)
  { return (6.0*t - 4.0); }
  inline double d2p2 (double t)
  { return (-d2p1(t)); }
  inline double d2q2 (double t)
  { return (6.0*t - 2.0); } 

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int dest, int i);
  void UpdateY (int source, int dest, int i);
  void UpdateZ (int source, int dest, int i);
  void UpdateXPeriodic (int source, int dest, int i);
  void UpdateYPeriodic (int source, int dest, int i);
  void UpdateZPeriodic (int source, int dest, int i);
  bool UpToDate, Periodic;
public:
  Array<TinyVector<double,16>,4> F;

  int Nx, Ny, Nz, N;
  Grid *Xgrid, *Ygrid, *Zgrid;
  TinyVector<Grid*,3> Grids;
  void Update();
  inline complex<double> operator()(int ix, int iy, int iz, int i) const
  { return complex<double> (F(ix,iy,ix, i)[0], F(ix,iy,iz,i)[8]); }
  inline void Set (int ix, int iy, int iz, int i, complex<double> val)
    { UpToDate=false; F(ix, iy, iz, i)[0] = val.real(); F(ix,iy,iz,i)[8] = val.imag(); }
  inline void operator()(double x, double y, double z, 
			 Array<complex<double>,1> &vals);
  inline void d_dx      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dy      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dz      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dx2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dy2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dz2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdy   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dydz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void Grad      (double x, double y, double z, 
			 Array<cVec3,1> &grads);
  inline void ValGrad   (double x, double y, double z, 
			 Array<complex<double>,1> &vals, 
			 Array<cVec3,  1> &grads);
  inline void Laplacian (double x, double y, double z, 
			 Array<complex<double>,1> &vals);

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid, int n)
  {
    Xgrid = xgrid; Nx = xgrid->NumPoints;
    Ygrid = ygrid; Ny = ygrid->NumPoints;
    Zgrid = zgrid; Nz = zgrid->NumPoints;
    N = n;
    
    F.resize(Nx,Ny,Nz,N);
    UpToDate = false;
  }
  
  /// Copy constructor
  inline ComplexMultiTricubicSpline (const ComplexMultiTricubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline ComplexMultiTricubicSpline & operator= (ComplexMultiTricubicSpline &a);
  inline ComplexMultiTricubicSpline & operator= (ComplexMultiTricubicSpline a);


  inline void Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
		    const Array<complex<double>,4> &init, bool periodic=false);

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid,
		      const Array<complex<double>,4> &init, bool periodic=false)
  {
    Init (xgrid, ygrid, zgrid, init, periodic);
  }
  ComplexMultiTricubicSpline() : UpToDate(false) 
  { /* Do nothing. */ }
};

inline 
ComplexMultiTricubicSpline::ComplexMultiTricubicSpline
(const ComplexMultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
}


inline 
ComplexMultiTricubicSpline& ComplexMultiTricubicSpline::operator=
(ComplexMultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}

inline ComplexMultiTricubicSpline& 
ComplexMultiTricubicSpline::operator=(ComplexMultiTricubicSpline a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}



inline void 
ComplexMultiTricubicSpline::Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
				  const Array<complex<double>,4> &init,
				  bool periodic)
{
  Periodic = periodic;

  Xgrid = xgrid; Nx = xgrid->NumPoints;
  Ygrid = ygrid; Ny = ygrid->NumPoints;
  Zgrid = zgrid; Nz = zgrid->NumPoints;

  assert (init.extent(0) == Nx);
  assert (init.extent(1) == Ny);
  assert (init.extent(2) == Nz);
  N = init.extent(3);

  F.resize(Nx,Ny,Nz,N);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int i=0; i<N; i++) {
	  F(ix,iy,iz,i)[0] = init(ix,iy,iz,i).real();
	  F(ix,iy,iz,i)[8] = init(ix,iy,iz,i).imag();
	}
  UpToDate = false;
  Update();
}




inline void ComplexMultiTricubicSpline::operator()
  (double x, double y, double z, Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);


  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;
  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);

  for (int i=0; i<N; i++) {
    //////////////////
    /// Real parts ///
    //////////////////
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double re, im;
    re = 
      a0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im = 
      a0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
    vals(i) = complex<double>(re,im);
  }
}


inline void 
ComplexMultiTricubicSpline::d_dx (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double re, im;
    re = 
      da0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      da1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      da2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      da3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));

    im = 
      da0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      da1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      da2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      da3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
    vals(i) = complex<double>(re,im);
  }
}



inline void
ComplexMultiTricubicSpline::d_dy (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);

  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz

    double re, im;
    re =
      a0*
      (db0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       db1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       db2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       db3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (db0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       db1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       db2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       db3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (db0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       db1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       db2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       db3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (db0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       db1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       db2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       db3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im =
      a0*
      (db0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       db1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       db2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       db3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (db0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       db1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       db2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       db3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (db0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       db1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       db2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       db3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (db0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       db1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       db2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       db3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
    vals(i) = complex<double> (re, im);
  }
}

inline void 
ComplexMultiTricubicSpline::d_dz (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  

  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz

    double re, im;
    re =
      a0*
      (b0*(rY000*dc0+rY001*dc1+rY002*dc2+rY003*dc3) +
       b1*(rY010*dc0+rY011*dc1+rY012*dc2+rY013*dc3) +
       b2*(rY020*dc0+rY021*dc1+rY022*dc2+rY023*dc3) +
       b3*(rY030*dc0+rY031*dc1+rY032*dc2+rY033*dc3))+
      a1 *
      (b0*(rY100*dc0+rY101*dc1+rY102*dc2+rY103*dc3) +
       b1*(rY110*dc0+rY111*dc1+rY112*dc2+rY113*dc3) +
       b2*(rY120*dc0+rY121*dc1+rY122*dc2+rY123*dc3) +
       b3*(rY130*dc0+rY131*dc1+rY132*dc2+rY133*dc3))+
      a2 *
      (b0*(rY200*dc0+rY201*dc1+rY202*dc2+rY203*dc3) +
       b1*(rY210*dc0+rY211*dc1+rY212*dc2+rY213*dc3) +
       b2*(rY220*dc0+rY221*dc1+rY222*dc2+rY223*dc3) +
       b3*(rY230*dc0+rY231*dc1+rY232*dc2+rY233*dc3))+
      a3 *
      (b0*(rY300*dc0+rY301*dc1+rY302*dc2+rY303*dc3) +
       b1*(rY310*dc0+rY311*dc1+rY312*dc2+rY313*dc3) +
       b2*(rY320*dc0+rY321*dc1+rY322*dc2+rY323*dc3) +
       b3*(rY330*dc0+rY331*dc1+rY332*dc2+rY333*dc3));
    im = 
      a0*
      (b0*(iY000*dc0+iY001*dc1+iY002*dc2+iY003*dc3) +
       b1*(iY010*dc0+iY011*dc1+iY012*dc2+iY013*dc3) +
       b2*(iY020*dc0+iY021*dc1+iY022*dc2+iY023*dc3) +
       b3*(iY030*dc0+iY031*dc1+iY032*dc2+iY033*dc3))+
      a1 *
      (b0*(iY100*dc0+iY101*dc1+iY102*dc2+iY103*dc3) +
       b1*(iY110*dc0+iY111*dc1+iY112*dc2+iY113*dc3) +
       b2*(iY120*dc0+iY121*dc1+iY122*dc2+iY123*dc3) +
       b3*(iY130*dc0+iY131*dc1+iY132*dc2+iY133*dc3))+
      a2 *
      (b0*(iY200*dc0+iY201*dc1+iY202*dc2+iY203*dc3) +
       b1*(iY210*dc0+iY211*dc1+iY212*dc2+iY213*dc3) +
       b2*(iY220*dc0+iY221*dc1+iY222*dc2+iY223*dc3) +
       b3*(iY230*dc0+iY231*dc1+iY232*dc2+iY233*dc3))+
      a3 *
      (b0*(iY300*dc0+iY301*dc1+iY302*dc2+iY303*dc3) +
       b1*(iY310*dc0+iY311*dc1+iY312*dc2+iY313*dc3) +
       b2*(iY320*dc0+iY321*dc1+iY322*dc2+iY323*dc3) +
       b3*(iY330*dc0+iY331*dc1+iY332*dc2+iY333*dc3));
    vals(i) = complex<double> (re, im);
  }
}

inline void
ComplexMultiTricubicSpline::Grad (double x, double y, double z, 
				  Array<cVec3,1> &grads)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz

    TinyVector<double,3> reGrad, imGrad;
    reGrad[0] = 
      da0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      da1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      da2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      da3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    imGrad[0] = 
      da0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      da1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      da2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      da3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    reGrad[1] = 
      a0*
      (db0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       db1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       db2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       db3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (db0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       db1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       db2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       db3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (db0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       db1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       db2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       db3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (db0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       db1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       db2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       db3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    imGrad[1] = 
      a0*
      (db0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       db1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       db2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       db3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (db0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       db1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       db2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       db3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (db0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       db1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       db2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       db3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (db0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       db1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       db2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       db3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    reGrad[2] = 
      a0*
      (b0*(rY000*dc0+rY001*dc1+rY002*dc2+rY003*dc3) +
       b1*(rY010*dc0+rY011*dc1+rY012*dc2+rY013*dc3) +
       b2*(rY020*dc0+rY021*dc1+rY022*dc2+rY023*dc3) +
       b3*(rY030*dc0+rY031*dc1+rY032*dc2+rY033*dc3))+
      a1 *
      (b0*(rY100*dc0+rY101*dc1+rY102*dc2+rY103*dc3) +
       b1*(rY110*dc0+rY111*dc1+rY112*dc2+rY113*dc3) +
       b2*(rY120*dc0+rY121*dc1+rY122*dc2+rY123*dc3) +
       b3*(rY130*dc0+rY131*dc1+rY132*dc2+rY133*dc3))+
      a2 *
      (b0*(rY200*dc0+rY201*dc1+rY202*dc2+rY203*dc3) +
       b1*(rY210*dc0+rY211*dc1+rY212*dc2+rY213*dc3) +
       b2*(rY220*dc0+rY221*dc1+rY222*dc2+rY223*dc3) +
       b3*(rY230*dc0+rY231*dc1+rY232*dc2+rY233*dc3))+
      a3 *
      (b0*(rY300*dc0+rY301*dc1+rY302*dc2+rY303*dc3) +
       b1*(rY310*dc0+rY311*dc1+rY312*dc2+rY313*dc3) +
       b2*(rY320*dc0+rY321*dc1+rY322*dc2+rY323*dc3) +
       b3*(rY330*dc0+rY331*dc1+rY332*dc2+rY333*dc3));
    imGrad[2] = 
      a0*
      (b0*(iY000*dc0+iY001*dc1+iY002*dc2+iY003*dc3) +
       b1*(iY010*dc0+iY011*dc1+iY012*dc2+iY013*dc3) +
       b2*(iY020*dc0+iY021*dc1+iY022*dc2+iY023*dc3) +
       b3*(iY030*dc0+iY031*dc1+iY032*dc2+iY033*dc3))+
      a1 *
      (b0*(iY100*dc0+iY101*dc1+iY102*dc2+iY103*dc3) +
       b1*(iY110*dc0+iY111*dc1+iY112*dc2+iY113*dc3) +
       b2*(iY120*dc0+iY121*dc1+iY122*dc2+iY123*dc3) +
       b3*(iY130*dc0+iY131*dc1+iY132*dc2+iY133*dc3))+
      a2 *
      (b0*(iY200*dc0+iY201*dc1+iY202*dc2+iY203*dc3) +
       b1*(iY210*dc0+iY211*dc1+iY212*dc2+iY213*dc3) +
       b2*(iY220*dc0+iY221*dc1+iY222*dc2+iY223*dc3) +
       b3*(iY230*dc0+iY231*dc1+iY232*dc2+iY233*dc3))+
      a3 *
      (b0*(iY300*dc0+iY301*dc1+iY302*dc2+iY303*dc3) +
       b1*(iY310*dc0+iY311*dc1+iY312*dc2+iY313*dc3) +
       b2*(iY320*dc0+iY321*dc1+iY322*dc2+iY323*dc3) +
       b3*(iY330*dc0+iY331*dc1+iY332*dc2+iY333*dc3));
    grads(i)[0] = complex<double> (reGrad[0], imGrad[0]);
    grads(i)[1] = complex<double> (reGrad[1], imGrad[1]);
    grads(i)[2] = complex<double> (reGrad[2], imGrad[2]);
  }
}

/// Returns the value and computes the gradient
inline void
ComplexMultiTricubicSpline::ValGrad(double x, double y, double z, 
				    Array<complex<double>,1> &vals, 
				    Array<cVec3,1> &grads)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  
    double re, im;
    Vec3 reGrad, imGrad;

    re = 
      a0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));

    im = 
      a0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    vals(i) = complex<double> (re, im);


    reGrad[0] = 
      da0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      da1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      da2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      da3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    imGrad[0] = 
      da0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      da1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      da2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      da3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    reGrad[1] = 
      a0*
      (db0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       db1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       db2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       db3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (db0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       db1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       db2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       db3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (db0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       db1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       db2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       db3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (db0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       db1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       db2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       db3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    imGrad[1] = 
      a0*
      (db0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       db1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       db2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       db3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (db0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       db1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       db2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       db3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (db0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       db1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       db2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       db3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (db0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       db1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       db2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       db3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    reGrad[2] = 
      a0*
      (b0*(rY000*dc0+rY001*dc1+rY002*dc2+rY003*dc3) +
       b1*(rY010*dc0+rY011*dc1+rY012*dc2+rY013*dc3) +
       b2*(rY020*dc0+rY021*dc1+rY022*dc2+rY023*dc3) +
       b3*(rY030*dc0+rY031*dc1+rY032*dc2+rY033*dc3))+
      a1 *
      (b0*(rY100*dc0+rY101*dc1+rY102*dc2+rY103*dc3) +
       b1*(rY110*dc0+rY111*dc1+rY112*dc2+rY113*dc3) +
       b2*(rY120*dc0+rY121*dc1+rY122*dc2+rY123*dc3) +
       b3*(rY130*dc0+rY131*dc1+rY132*dc2+rY133*dc3))+
      a2 *
      (b0*(rY200*dc0+rY201*dc1+rY202*dc2+rY203*dc3) +
       b1*(rY210*dc0+rY211*dc1+rY212*dc2+rY213*dc3) +
       b2*(rY220*dc0+rY221*dc1+rY222*dc2+rY223*dc3) +
       b3*(rY230*dc0+rY231*dc1+rY232*dc2+rY233*dc3))+
      a3 *
      (b0*(rY300*dc0+rY301*dc1+rY302*dc2+rY303*dc3) +
       b1*(rY310*dc0+rY311*dc1+rY312*dc2+rY313*dc3) +
       b2*(rY320*dc0+rY321*dc1+rY322*dc2+rY323*dc3) +
       b3*(rY330*dc0+rY331*dc1+rY332*dc2+rY333*dc3));
    imGrad[2] = 
      a0*
      (b0*(iY000*dc0+iY001*dc1+iY002*dc2+iY003*dc3) +
       b1*(iY010*dc0+iY011*dc1+iY012*dc2+iY013*dc3) +
       b2*(iY020*dc0+iY021*dc1+iY022*dc2+iY023*dc3) +
       b3*(iY030*dc0+iY031*dc1+iY032*dc2+iY033*dc3))+
      a1 *
      (b0*(iY100*dc0+iY101*dc1+iY102*dc2+iY103*dc3) +
       b1*(iY110*dc0+iY111*dc1+iY112*dc2+iY113*dc3) +
       b2*(iY120*dc0+iY121*dc1+iY122*dc2+iY123*dc3) +
       b3*(iY130*dc0+iY131*dc1+iY132*dc2+iY133*dc3))+
      a2 *
      (b0*(iY200*dc0+iY201*dc1+iY202*dc2+iY203*dc3) +
       b1*(iY210*dc0+iY211*dc1+iY212*dc2+iY213*dc3) +
       b2*(iY220*dc0+iY221*dc1+iY222*dc2+iY223*dc3) +
       b3*(iY230*dc0+iY231*dc1+iY232*dc2+iY233*dc3))+
      a3 *
      (b0*(iY300*dc0+iY301*dc1+iY302*dc2+iY303*dc3) +
       b1*(iY310*dc0+iY311*dc1+iY312*dc2+iY313*dc3) +
       b2*(iY320*dc0+iY321*dc1+iY322*dc2+iY323*dc3) +
       b3*(iY330*dc0+iY331*dc1+iY332*dc2+iY333*dc3));
    grads(i)[0] = complex<double>(reGrad[0], imGrad[0]);
    grads(i)[1] = complex<double>(reGrad[1], imGrad[1]);
    grads(i)[2] = complex<double>(reGrad[2], imGrad[2]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dx2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  
    double re, im;

    re = 
      d2a0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      d2a1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      d2a2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      d2a3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));

    im = 
      d2a0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      d2a1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      d2a2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      d2a3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    vals(i) = complex<double> (re, im);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dy2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  
    double re, im;
    re = 
      a0*
      (d2b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       d2b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       d2b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       d2b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (d2b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       d2b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       d2b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       d2b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (d2b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       d2b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       d2b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       d2b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (d2b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       d2b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       d2b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       d2b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im = 
      a0*
      (d2b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       d2b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       d2b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       d2b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (d2b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       d2b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       d2b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       d2b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (d2b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       d2b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       d2b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       d2b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (d2b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       d2b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       d2b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       d2b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
      vals(i) = complex<double> (re, im);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dz2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  
    double re, im;
    re = 
      a0*
      (b0*(rY000*d2c0+rY001*d2c1+rY002*d2c2+rY003*d2c3) +
       b1*(rY010*d2c0+rY011*d2c1+rY012*d2c2+rY013*d2c3) +
       b2*(rY020*d2c0+rY021*d2c1+rY022*d2c2+rY023*d2c3) +
       b3*(rY030*d2c0+rY031*d2c1+rY032*d2c2+rY033*d2c3))+
      a1 *
      (b0*(rY100*d2c0+rY101*d2c1+rY102*d2c2+rY103*d2c3) +
       b1*(rY110*d2c0+rY111*d2c1+rY112*d2c2+rY113*d2c3) +
       b2*(rY120*d2c0+rY121*d2c1+rY122*d2c2+rY123*d2c3) +
       b3*(rY130*d2c0+rY131*d2c1+rY132*d2c2+rY133*d2c3))+
      a2 *
      (b0*(rY200*d2c0+rY201*d2c1+rY202*d2c2+rY203*d2c3) +
       b1*(rY210*d2c0+rY211*d2c1+rY212*d2c2+rY213*d2c3) +
       b2*(rY220*d2c0+rY221*d2c1+rY222*d2c2+rY223*d2c3) +
       b3*(rY230*d2c0+rY231*d2c1+rY232*d2c2+rY233*d2c3))+
      a3 *
      (b0*(rY300*d2c0+rY301*d2c1+rY302*d2c2+rY303*d2c3) +
       b1*(rY310*d2c0+rY311*d2c1+rY312*d2c2+rY313*d2c3) +
       b2*(rY320*d2c0+rY321*d2c1+rY322*d2c2+rY323*d2c3) +
       b3*(rY330*d2c0+rY331*d2c1+rY332*d2c2+rY333*d2c3));
    im = 
      a0*
      (b0*(iY000*d2c0+iY001*d2c1+iY002*d2c2+iY003*d2c3) +
       b1*(iY010*d2c0+iY011*d2c1+iY012*d2c2+iY013*d2c3) +
       b2*(iY020*d2c0+iY021*d2c1+iY022*d2c2+iY023*d2c3) +
       b3*(iY030*d2c0+iY031*d2c1+iY032*d2c2+iY033*d2c3))+
      a1 *
      (b0*(iY100*d2c0+iY101*d2c1+iY102*d2c2+iY103*d2c3) +
       b1*(iY110*d2c0+iY111*d2c1+iY112*d2c2+iY113*d2c3) +
       b2*(iY120*d2c0+iY121*d2c1+iY122*d2c2+iY123*d2c3) +
       b3*(iY130*d2c0+iY131*d2c1+iY132*d2c2+iY133*d2c3))+
      a2 *
      (b0*(iY200*d2c0+iY201*d2c1+iY202*d2c2+iY203*d2c3) +
       b1*(iY210*d2c0+iY211*d2c1+iY212*d2c2+iY213*d2c3) +
       b2*(iY220*d2c0+iY221*d2c1+iY222*d2c2+iY223*d2c3) +
       b3*(iY230*d2c0+iY231*d2c1+iY232*d2c2+iY233*d2c3))+
      a3 *
      (b0*(iY300*d2c0+iY301*d2c1+iY302*d2c2+iY303*d2c3) +
       b1*(iY310*d2c0+iY311*d2c1+iY312*d2c2+iY313*d2c3) +
       b2*(iY320*d2c0+iY321*d2c1+iY322*d2c2+iY323*d2c3) +
       b3*(iY330*d2c0+iY331*d2c1+iY332*d2c2+iY333*d2c3));
    vals(i) = complex<double>(re, im);
  }
}


inline void 
ComplexMultiTricubicSpline::Laplacian (double x, double y, double z, 
				       Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  
    double re, im;

    re = 
      d2a0*
      (b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      d2a1 *
      (b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      d2a2 *
      (b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      d2a3 *
      (b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im = 
      d2a0*
      (b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      d2a1 *
      (b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      d2a2 *
      (b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      d2a3 *
      (b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
  
    re +=
      a0*
      (d2b0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       d2b1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       d2b2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       d2b3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      a1 *
      (d2b0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       d2b1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       d2b2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       d2b3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      a2 *
      (d2b0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       d2b1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       d2b2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       d2b3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      a3 *
      (d2b0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       d2b1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       d2b2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       d2b3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im +=
      a0*
      (d2b0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       d2b1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       d2b2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       d2b3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      a1 *
      (d2b0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       d2b1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       d2b2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       d2b3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      a2 *
      (d2b0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       d2b1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       d2b2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       d2b3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      a3 *
      (d2b0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       d2b1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       d2b2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       d2b3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));

    re += 
      a0*
      (b0*(rY000*d2c0+rY001*d2c1+rY002*d2c2+rY003*d2c3) +
       b1*(rY010*d2c0+rY011*d2c1+rY012*d2c2+rY013*d2c3) +
       b2*(rY020*d2c0+rY021*d2c1+rY022*d2c2+rY023*d2c3) +
       b3*(rY030*d2c0+rY031*d2c1+rY032*d2c2+rY033*d2c3))+
      a1 *
      (b0*(rY100*d2c0+rY101*d2c1+rY102*d2c2+rY103*d2c3) +
       b1*(rY110*d2c0+rY111*d2c1+rY112*d2c2+rY113*d2c3) +
       b2*(rY120*d2c0+rY121*d2c1+rY122*d2c2+rY123*d2c3) +
       b3*(rY130*d2c0+rY131*d2c1+rY132*d2c2+rY133*d2c3))+
      a2 *
      (b0*(rY200*d2c0+rY201*d2c1+rY202*d2c2+rY203*d2c3) +
       b1*(rY210*d2c0+rY211*d2c1+rY212*d2c2+rY213*d2c3) +
       b2*(rY220*d2c0+rY221*d2c1+rY222*d2c2+rY223*d2c3) +
       b3*(rY230*d2c0+rY231*d2c1+rY232*d2c2+rY233*d2c3))+
      a3 *
      (b0*(rY300*d2c0+rY301*d2c1+rY302*d2c2+rY303*d2c3) +
       b1*(rY310*d2c0+rY311*d2c1+rY312*d2c2+rY313*d2c3) +
       b2*(rY320*d2c0+rY321*d2c1+rY322*d2c2+rY323*d2c3) +
       b3*(rY330*d2c0+rY331*d2c1+rY332*d2c2+rY333*d2c3));
    im += 
      a0*
      (b0*(iY000*d2c0+iY001*d2c1+iY002*d2c2+iY003*d2c3) +
       b1*(iY010*d2c0+iY011*d2c1+iY012*d2c2+iY013*d2c3) +
       b2*(iY020*d2c0+iY021*d2c1+iY022*d2c2+iY023*d2c3) +
       b3*(iY030*d2c0+iY031*d2c1+iY032*d2c2+iY033*d2c3))+
      a1 *
      (b0*(iY100*d2c0+iY101*d2c1+iY102*d2c2+iY103*d2c3) +
       b1*(iY110*d2c0+iY111*d2c1+iY112*d2c2+iY113*d2c3) +
       b2*(iY120*d2c0+iY121*d2c1+iY122*d2c2+iY123*d2c3) +
       b3*(iY130*d2c0+iY131*d2c1+iY132*d2c2+iY133*d2c3))+
      a2 *
      (b0*(iY200*d2c0+iY201*d2c1+iY202*d2c2+iY203*d2c3) +
       b1*(iY210*d2c0+iY211*d2c1+iY212*d2c2+iY213*d2c3) +
       b2*(iY220*d2c0+iY221*d2c1+iY222*d2c2+iY223*d2c3) +
       b3*(iY230*d2c0+iY231*d2c1+iY232*d2c2+iY233*d2c3))+
      a3 *
      (b0*(iY300*d2c0+iY301*d2c1+iY302*d2c2+iY303*d2c3) +
       b1*(iY310*d2c0+iY311*d2c1+iY312*d2c2+iY313*d2c3) +
       b2*(iY320*d2c0+iY321*d2c1+iY322*d2c2+iY323*d2c3) +
       b3*(iY330*d2c0+iY331*d2c1+iY332*d2c2+iY333*d2c3));
    vals(i) = complex<double>(re,im);
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdy (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double re, im;
    re =
      da0*
      (db0*(rY000*c0+rY001*c1+rY002*c2+rY003*c3) +
       db1*(rY010*c0+rY011*c1+rY012*c2+rY013*c3) +
       db2*(rY020*c0+rY021*c1+rY022*c2+rY023*c3) +
       db3*(rY030*c0+rY031*c1+rY032*c2+rY033*c3))+
      da1 *
      (db0*(rY100*c0+rY101*c1+rY102*c2+rY103*c3) +
       db1*(rY110*c0+rY111*c1+rY112*c2+rY113*c3) +
       db2*(rY120*c0+rY121*c1+rY122*c2+rY123*c3) +
       db3*(rY130*c0+rY131*c1+rY132*c2+rY133*c3))+
      da2 *
      (db0*(rY200*c0+rY201*c1+rY202*c2+rY203*c3) +
       db1*(rY210*c0+rY211*c1+rY212*c2+rY213*c3) +
       db2*(rY220*c0+rY221*c1+rY222*c2+rY223*c3) +
       db3*(rY230*c0+rY231*c1+rY232*c2+rY233*c3))+
      da3 *
      (db0*(rY300*c0+rY301*c1+rY302*c2+rY303*c3) +
       db1*(rY310*c0+rY311*c1+rY312*c2+rY313*c3) +
       db2*(rY320*c0+rY321*c1+rY322*c2+rY323*c3) +
       db3*(rY330*c0+rY331*c1+rY332*c2+rY333*c3));
    im =
      da0*
      (db0*(iY000*c0+iY001*c1+iY002*c2+iY003*c3) +
       db1*(iY010*c0+iY011*c1+iY012*c2+iY013*c3) +
       db2*(iY020*c0+iY021*c1+iY022*c2+iY023*c3) +
       db3*(iY030*c0+iY031*c1+iY032*c2+iY033*c3))+
      da1 *
      (db0*(iY100*c0+iY101*c1+iY102*c2+iY103*c3) +
       db1*(iY110*c0+iY111*c1+iY112*c2+iY113*c3) +
       db2*(iY120*c0+iY121*c1+iY122*c2+iY123*c3) +
       db3*(iY130*c0+iY131*c1+iY132*c2+iY133*c3))+
      da2 *
      (db0*(iY200*c0+iY201*c1+iY202*c2+iY203*c3) +
       db1*(iY210*c0+iY211*c1+iY212*c2+iY213*c3) +
       db2*(iY220*c0+iY221*c1+iY222*c2+iY223*c3) +
       db3*(iY230*c0+iY231*c1+iY232*c2+iY233*c3))+
      da3 *
      (db0*(iY300*c0+iY301*c1+iY302*c2+iY303*c3) +
       db1*(iY310*c0+iY311*c1+iY312*c2+iY313*c3) +
       db2*(iY320*c0+iY321*c1+iY322*c2+iY323*c3) +
       db3*(iY330*c0+iY331*c1+iY332*c2+iY333*c3));
    vals(i) = complex<double>(re,im);
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);


  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz
  

    double re, im;
    re = 
      da0*
      (b0*(rY000*dc0+rY001*dc1+rY002*dc2+rY003*dc3) +
       b1*(rY010*dc0+rY011*dc1+rY012*dc2+rY013*dc3) +
       b2*(rY020*dc0+rY021*dc1+rY022*dc2+rY023*dc3) +
       b3*(rY030*dc0+rY031*dc1+rY032*dc2+rY033*dc3))+
      da1 *
      (b0*(rY100*dc0+rY101*dc1+rY102*dc2+rY103*dc3) +
       b1*(rY110*dc0+rY111*dc1+rY112*dc2+rY113*dc3) +
       b2*(rY120*dc0+rY121*dc1+rY122*dc2+rY123*dc3) +
       b3*(rY130*dc0+rY131*dc1+rY132*dc2+rY133*dc3))+
      da2 *
      (b0*(rY200*dc0+rY201*dc1+rY202*dc2+rY203*dc3) +
       b1*(rY210*dc0+rY211*dc1+rY212*dc2+rY213*dc3) +
       b2*(rY220*dc0+rY221*dc1+rY222*dc2+rY223*dc3) +
       b3*(rY230*dc0+rY231*dc1+rY232*dc2+rY233*dc3))+
      da3 *
      (b0*(rY300*dc0+rY301*dc1+rY302*dc2+rY303*dc3) +
       b1*(rY310*dc0+rY311*dc1+rY312*dc2+rY313*dc3) +
       b2*(rY320*dc0+rY321*dc1+rY322*dc2+rY323*dc3) +
       b3*(rY330*dc0+rY331*dc1+rY332*dc2+rY333*dc3));
    im = 
      da0*
      (b0*(iY000*dc0+iY001*dc1+iY002*dc2+iY003*dc3) +
       b1*(iY010*dc0+iY011*dc1+iY012*dc2+iY013*dc3) +
       b2*(iY020*dc0+iY021*dc1+iY022*dc2+iY023*dc3) +
       b3*(iY030*dc0+iY031*dc1+iY032*dc2+iY033*dc3))+
      da1 *
      (b0*(iY100*dc0+iY101*dc1+iY102*dc2+iY103*dc3) +
       b1*(iY110*dc0+iY111*dc1+iY112*dc2+iY113*dc3) +
       b2*(iY120*dc0+iY121*dc1+iY122*dc2+iY123*dc3) +
       b3*(iY130*dc0+iY131*dc1+iY132*dc2+iY133*dc3))+
      da2 *
      (b0*(iY200*dc0+iY201*dc1+iY202*dc2+iY203*dc3) +
       b1*(iY210*dc0+iY211*dc1+iY212*dc2+iY213*dc3) +
       b2*(iY220*dc0+iY221*dc1+iY222*dc2+iY223*dc3) +
       b3*(iY230*dc0+iY231*dc1+iY232*dc2+iY233*dc3))+
      da3 *
      (b0*(iY300*dc0+iY301*dc1+iY302*dc2+iY303*dc3) +
       b1*(iY310*dc0+iY311*dc1+iY312*dc2+iY313*dc3) +
       b2*(iY320*dc0+iY321*dc1+iY322*dc2+iY323*dc3) +
       b3*(iY330*dc0+iY331*dc1+iY332*dc2+iY333*dc3));
    vals(i) = complex<double>(re,im);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dydz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);

  for (int i=0; i<N; i++) {
    double& rY000 = F(ix,iy,iz,i)[0];      //   F
    double& rY001 = F(ix,iy,iz+1,i)[0];    //   F
    double& rY002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& rY003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& rY010 = F(ix,iy+1,iz,i)[0];    //   F
    double& rY011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& rY012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& rY013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& rY021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& rY022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& rY023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& rY030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& rY031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& rY033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY100 = F(ix+1,iy,iz,i)[0];      //   F
    double& rY101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& rY102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& rY103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& rY110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& rY111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& rY112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& rY113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& rY120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& rY121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& rY122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& rY123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& rY130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& rY131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& rY132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& rY133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& rY200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& rY201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& rY202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& rY203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& rY211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& rY221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& rY223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& rY300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& rY301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& rY302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& rY303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& rY310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& rY311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& rY312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& rY313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& rY320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& rY321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& rY322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& rY323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& rY330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& rY331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& rY332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& rY333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    ///////////////////////
    /// Imaginary parts ///
    ///////////////////////
    double& iY000 = F(ix,iy,iz,i)[8];      //   F
    double& iY001 = F(ix,iy,iz+1,i)[8];    //   F
    double& iY002 = F(ix,iy,iz,i)[11];      //  dF/dz
    double& iY003 = F(ix,iy,iz+1,i)[11];    //  dF/dz
    double& iY010 = F(ix,iy+1,iz,i)[8];    //   F
    double& iY011 = F(ix,iy+1,iz+1,i)[8];  //   F
    double& iY012 = F(ix,iy+1,iz,i)[11];    //  dF/dz
    double& iY013 = F(ix,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY020 = F(ix,iy,iz,i)[10];      //  dF/dy
    double& iY021 = F(ix,iy,iz+1,i)[10];    //  dF/dy
    double& iY022 = F(ix,iy,iz,i)[14];      // d2F/dydz
    double& iY023 = F(ix,iy,iz+1,i)[14];    // d2F/dydz
    double& iY030 = F(ix,iy+1,iz,i)[10];    //  dF/dy
    double& iY031 = F(ix,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY032 = F(ix,iy+1,iz,i)[14];    // d2F/dydz
    double& iY033 = F(ix,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY100 = F(ix+1,iy,iz,i)[8];      //   F
    double& iY101 = F(ix+1,iy,iz+1,i)[8];    //   F
    double& iY102 = F(ix+1,iy,iz,i)[11];      //  dF/dz
    double& iY103 = F(ix+1,iy,iz+1,i)[11];    //  dF/dz
    double& iY110 = F(ix+1,iy+1,iz,i)[8];    //   F
    double& iY111 = F(ix+1,iy+1,iz+1,i)[8];  //   F
    double& iY112 = F(ix+1,iy+1,iz,i)[11];    //  dF/dz
    double& iY113 = F(ix+1,iy+1,iz+1,i)[11];  //  dF/dz
    double& iY120 = F(ix+1,iy,iz,i)[10];      //  dF/dy
    double& iY121 = F(ix+1,iy,iz+1,i)[10];    //  dF/dy
    double& iY122 = F(ix+1,iy,iz,i)[14];      // d2F/dydz
    double& iY123 = F(ix+1,iy,iz+1,i)[14];    // d2F/dydz
    double& iY130 = F(ix+1,iy+1,iz,i)[10];    //  dF/dy
    double& iY131 = F(ix+1,iy+1,iz+1,i)[10];  //  dF/dy
    double& iY132 = F(ix+1,iy+1,iz,i)[14];    // d2F/dydz
    double& iY133 = F(ix+1,iy+1,iz+1,i)[14];  // d2F/dydz
    
    double& iY200 = F(ix,iy,iz,i)[9];      //  dF/dx
    double& iY201 = F(ix,iy,iz+1,i)[9];    //  dF/dx
    double& iY202 = F(ix,iy,iz,i)[13];      // d2F/dxdz
    double& iY203 = F(ix,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY210 = F(ix,iy+1,iz,i)[9];    //  dF/dx
    double& iY211 = F(ix,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY212 = F(ix,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY213 = F(ix,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY220 = F(ix,iy,iz,i)[12];      // d2F/dxdy
    double& iY221 = F(ix,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY222 = F(ix,iy,iz,i)[15];      // d3F/dxdydz
    double& iY223 = F(ix,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY230 = F(ix,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY231 = F(ix,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY232 = F(ix,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY233 = F(ix,iy+1,iz+1,i)[15];  // d3F/dxdydz
    
    double& iY300 = F(ix+1,iy,iz,i)[9];      //  dF/dx
    double& iY301 = F(ix+1,iy,iz+1,i)[9];    //  dF/dx
    double& iY302 = F(ix+1,iy,iz,i)[13];      // d2F/dxdz
    double& iY303 = F(ix+1,iy,iz+1,i)[13];    // d2F/dxdz
    double& iY310 = F(ix+1,iy+1,iz,i)[9];    //  dF/dx
    double& iY311 = F(ix+1,iy+1,iz+1,i)[9];  //  dF/dx
    double& iY312 = F(ix+1,iy+1,iz,i)[13];    // d2F/dxdz
    double& iY313 = F(ix+1,iy+1,iz+1,i)[13];  // d2F/dxdz
    double& iY320 = F(ix+1,iy,iz,i)[12];      // d2F/dxdy
    double& iY321 = F(ix+1,iy,iz+1,i)[12];    // d2F/dxdy
    double& iY322 = F(ix+1,iy,iz,i)[15];      // d3F/dxdydz
    double& iY323 = F(ix+1,iy,iz+1,i)[15];    // d3F/dxdydz
    double& iY330 = F(ix+1,iy+1,iz,i)[12];    // d2F/dxdy
    double& iY331 = F(ix+1,iy+1,iz+1,i)[12];  // d2F/dxdy
    double& iY332 = F(ix+1,iy+1,iz,i)[15];    // d3F/dxdydz
    double& iY333 = F(ix+1,iy+1,iz+1,i)[15];  // d3F/dxdydz

    double re, im;
    re = 
      a0*
      (db0*(rY000*dc0+rY001*dc1+rY002*dc2+rY003*dc3) +
       db1*(rY010*dc0+rY011*dc1+rY012*dc2+rY013*dc3) +
       db2*(rY020*dc0+rY021*dc1+rY022*dc2+rY023*dc3) +
       db3*(rY030*dc0+rY031*dc1+rY032*dc2+rY033*dc3))+
      a1 *
      (db0*(rY100*dc0+rY101*dc1+rY102*dc2+rY103*dc3) +
       db1*(rY110*dc0+rY111*dc1+rY112*dc2+rY113*dc3) +
       db2*(rY120*dc0+rY121*dc1+rY122*dc2+rY123*dc3) +
       db3*(rY130*dc0+rY131*dc1+rY132*dc2+rY133*dc3))+
      a2 *
      (db0*(rY200*dc0+rY201*dc1+rY202*dc2+rY203*dc3) +
       db1*(rY210*dc0+rY211*dc1+rY212*dc2+rY213*dc3) +
       db2*(rY220*dc0+rY221*dc1+rY222*dc2+rY223*dc3) +
       db3*(rY230*dc0+rY231*dc1+rY232*dc2+rY233*dc3))+
      a3 *
      (db0*(rY300*dc0+rY301*dc1+rY302*dc2+rY303*dc3) +
       db1*(rY310*dc0+rY311*dc1+rY312*dc2+rY313*dc3) +
       db2*(rY320*dc0+rY321*dc1+rY322*dc2+rY323*dc3) +
       db3*(rY330*dc0+rY331*dc1+rY332*dc2+rY333*dc3));
    im = 
      a0*
      (db0*(iY000*dc0+iY001*dc1+iY002*dc2+iY003*dc3) +
       db1*(iY010*dc0+iY011*dc1+iY012*dc2+iY013*dc3) +
       db2*(iY020*dc0+iY021*dc1+iY022*dc2+iY023*dc3) +
       db3*(iY030*dc0+iY031*dc1+iY032*dc2+iY033*dc3))+
      a1 *
      (db0*(iY100*dc0+iY101*dc1+iY102*dc2+iY103*dc3) +
       db1*(iY110*dc0+iY111*dc1+iY112*dc2+iY113*dc3) +
       db2*(iY120*dc0+iY121*dc1+iY122*dc2+iY123*dc3) +
       db3*(iY130*dc0+iY131*dc1+iY132*dc2+iY133*dc3))+
      a2 *
      (db0*(iY200*dc0+iY201*dc1+iY202*dc2+iY203*dc3) +
       db1*(iY210*dc0+iY211*dc1+iY212*dc2+iY213*dc3) +
       db2*(iY220*dc0+iY221*dc1+iY222*dc2+iY223*dc3) +
       db3*(iY230*dc0+iY231*dc1+iY232*dc2+iY233*dc3))+
      a3 *
      (db0*(iY300*dc0+iY301*dc1+iY302*dc2+iY303*dc3) +
       db1*(iY310*dc0+iY311*dc1+iY312*dc2+iY313*dc3) +
       db2*(iY320*dc0+iY321*dc1+iY322*dc2+iY323*dc3) +
       db3*(iY330*dc0+iY331*dc1+iY332*dc2+iY333*dc3));
      vals(i) = complex<double>(re,im);

  }
}

/// This replicates the first points of the 3D array to the last
/// points, making the function periodic
void MakePeriodic(Array<complex<double>,4> &A);

#endif
