#ifndef TRICUBIC_NUB_SPLINE_H
#define TRICUBIC_NUB_SPLINE_H

#include "BsplineHelper.h"
#include "NUBsplineBasis.h"
#include "Grid.h"

template<typename T, typename XGridType=LinearGrid, 
	 typename YGridType=XGridType, typename ZGridType=YGridType>
class TricubicNUBspline
{
private:
  NUBsplineBasis<XGridType> XBasis;
  NUBsplineBasis<YGridType> YBasis;
  NUBsplineBasis<ZGridType> ZBasis;
  int Nx, Ny, Nz;
  // The starting and ending values for the uniform grids
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  // The box dimensions and their inverses
  double Lx, LxInv, Ly, LyInv, Lz, LzInv;
  // The control points
  Array<T,3> P;

  TinyVector<double,3> Periodic;  
public:
  inline void Init (XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<T,3> &data, 
		    BCType xbct=PERIODIC, BCType ybct=PERIODIC, BCType zbct=PERIODIC);
  inline T operator() (TinyVector<double,3> r) const;
  inline void Evaluate (TinyVector<double,3> r, T &val, 
			TinyVector<T,3> &grad) const;
  inline void Evaluate (TinyVector<double,3> r, T &val,
			TinyVector<T,3> &grad, T &laplacian) const;
  inline void Evaluate (TinyVector<double,3> r, T & val,
			TinyVector<T,3> &grad, 
			TinyMatrix<T,3,3> &secDerivs) const;
};




void Duplicate (TinyVector<double,4> source,
		TinyVector<complex<double>,4> dest)
{
  dest[0] = complex<double>(source[0], source[0]);
  dest[1] = complex<double>(source[1], source[1]);
  dest[2] = complex<double>(source[2], source[2]);
  dest[3] = complex<double>(source[3], source[3]);
}

void Duplicate (TinyVector<double,4> source,
		TinyVector<double,4> dest)
{
  dest = source;
}



template<typename T, typename XGridType, typename YGridType, typename ZGridType>
void 
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Init
(XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<T,3> &data, 
 BCType xbc, BCType ybc, BCType zbc)
{
  // Set up 1D basis functions
  XBasis.Init (xgrid, xbc==PERIODIC);
  YBasis.Init (ygrid, ybc==PERIODIC);
  ZBasis.Init (zgrid, zbc==PERIODIC);
  Periodic[0] = xbc==PERIODIC ? 1.0 : 0.0;
  Periodic[1] = ybc==PERIODIC ? 1.0 : 0.0;
  Periodic[2] = zbc==PERIODIC ? 1.0 : 0.0;

  Nx = data.extent(0);
  Ny = data.extent(1);
  Nz = data.extent(2);

  int Mx, My, Mz;
  Mx = (xbc == PERIODIC) ? Nx+3 : Nx+2;
  My = (ybc == PERIODIC) ? Ny+3 : Ny+2;
  Mz = (zbc == PERIODIC) ? Nz+3 : Nz+2;
  assert (xgrid->NumPoints == Mx-2);
  assert (ygrid->NumPoints == My-2);
  assert (zgrid->NumPoints == Mz-2);

  P.resize(Mx, My, Mz);

  
  // Now solve interpolating equations
  TinyVector<T,4> lBC, rBC, dummy1, dummy2;
  ////////////////////
  // Do X direction //
  ////////////////////
  if (xbc == PERIODIC) {
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	SolvePeriodicInterp1D(XBasis, data(Range::all(), iy, iz), P(Range::all(),iy+1, iz+1));
  }
  else {
    if (xbc == FLAT) {
      XBasis (0   , dummy1, lBC);
      XBasis (Nx-1, dummy1, rBC);
    }
    else if (xbc == NATURAL) {
      XBasis (0   , dummy1, dummy2, lBC);
      XBasis (Nx-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	SolveDerivInterp1D(XBasis, data(Range::all(), iy, iz),  P(Range::all(), iy+1, iz+1), lBC, rBC);
  }
  ////////////////////
  // Do Y direction //
  ////////////////////
  if (ybc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++) 
	SolvePeriodicInterp1D(YBasis, P(ix, Range(1,Ny), iz), P(ix, Range::all(), iz));
  }
  else {
    if (ybc == FLAT) {
      YBasis (0,    dummy1, lBC);
      YBasis (Ny-1, dummy1, rBC);
    }
    else if (ybc == NATURAL) {
      YBasis (0,    dummy1, dummy2, lBC);
      YBasis (Ny-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++)
	SolveDerivInterp1D(YBasis, P(ix,Range(1,Ny), iz),  P(ix, Range::all(), iz), lBC, rBC);
  }
  ////////////////////
  // Do Z direction //
  ////////////////////
  if (zbc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++) 
	SolvePeriodicInterp1D(ZBasis, P(ix, iy, Range(1,Nz)), P(ix, iy, Range::all()));
  }
  else {
    if (zbc == FLAT) {
      ZBasis (0,    dummy1, lBC);
      ZBasis (Nz-1, dummy1, rBC);
    }
    else if (zbc == NATURAL) {
      ZBasis (0,    dummy1, dummy2, lBC);
      ZBasis (Nz-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++)
	SolveDerivInterp1D(ZBasis, P(ix, iy, Range(1,Nz)),  P(ix, iy, Range::all()), lBC, rBC);
  }
}

template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline T 
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::operator() (TinyVector<double,3> r) const
{
  TinyVector<double,4> a, b, c;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  // Return tensor product
  return 
    (a[0]*(b[0]*(c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3))+
	   b[1]*(c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3))+
	   b[2]*(c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3))+
	   b[3]*(c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3)))+
     a[1]*(b[0]*(c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3))+
	   b[1]*(c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3))+
	   b[2]*(c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3))+
	   b[3]*(c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3)))+
     a[2]*(b[0]*(c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3))+
	   b[1]*(c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3))+
	   b[2]*(c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3))+
	   b[3]*(c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3)))+
     a[3]*(b[0]*(c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3))+
	   b[1]*(c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3))+
	   b[2]*(c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3))+
	   b[3]*(c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3))));
}

template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Evaluate (TinyVector<double,3> r,
							      T &val,
							      TinyVector<T,3> &grad) const
{
  TinyVector<double,4> a, da, b, db, c, dc;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a, da); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b, db); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c, dc); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  T Pi[64], cP[16], bcP[4];
  Pi[ 0]=P(ix0,iy0,iz0);  Pi[ 1]=P(ix0,iy0,iz1);  Pi[ 2]=P(ix0,iy0,iz2);  Pi[ 3]=P(ix0,iy0,iz3);
  Pi[ 4]=P(ix0,iy1,iz0);  Pi[ 5]=P(ix0,iy1,iz1);  Pi[ 6]=P(ix0,iy1,iz2);  Pi[ 7]=P(ix0,iy1,iz3);
  Pi[ 8]=P(ix0,iy2,iz0);  Pi[ 9]=P(ix0,iy2,iz1);  Pi[10]=P(ix0,iy2,iz2);  Pi[11]=P(ix0,iy2,iz3);
  Pi[12]=P(ix0,iy3,iz0);  Pi[13]=P(ix0,iy3,iz1);  Pi[14]=P(ix0,iy3,iz2);  Pi[15]=P(ix0,iy3,iz3);
  Pi[16]=P(ix1,iy0,iz0);  Pi[17]=P(ix1,iy0,iz1);  Pi[18]=P(ix1,iy0,iz2);  Pi[19]=P(ix1,iy0,iz3);
  Pi[20]=P(ix1,iy1,iz0);  Pi[21]=P(ix1,iy1,iz1);  Pi[22]=P(ix1,iy1,iz2);  Pi[23]=P(ix1,iy1,iz3);
  Pi[24]=P(ix1,iy2,iz0);  Pi[25]=P(ix1,iy2,iz1);  Pi[26]=P(ix1,iy2,iz2);  Pi[27]=P(ix1,iy2,iz3);
  Pi[28]=P(ix1,iy3,iz0);  Pi[29]=P(ix1,iy3,iz1);  Pi[30]=P(ix1,iy3,iz2);  Pi[31]=P(ix1,iy3,iz3);
  Pi[32]=P(ix2,iy0,iz0);  Pi[33]=P(ix2,iy0,iz1);  Pi[34]=P(ix2,iy0,iz2);  Pi[35]=P(ix2,iy0,iz3);
  Pi[36]=P(ix2,iy1,iz0);  Pi[37]=P(ix2,iy1,iz1);  Pi[38]=P(ix2,iy1,iz2);  Pi[39]=P(ix2,iy1,iz3);
  Pi[40]=P(ix2,iy2,iz0);  Pi[41]=P(ix2,iy2,iz1);  Pi[42]=P(ix2,iy2,iz2);  Pi[43]=P(ix2,iy2,iz3);
  Pi[44]=P(ix2,iy3,iz0);  Pi[45]=P(ix2,iy3,iz1);  Pi[46]=P(ix2,iy3,iz2);  Pi[47]=P(ix2,iy3,iz3);
  Pi[48]=P(ix3,iy0,iz0);  Pi[49]=P(ix3,iy0,iz1);  Pi[50]=P(ix3,iy0,iz2);  Pi[51]=P(ix3,iy0,iz3);
  Pi[52]=P(ix3,iy1,iz0);  Pi[53]=P(ix3,iy1,iz1);  Pi[54]=P(ix3,iy1,iz2);  Pi[55]=P(ix3,iy1,iz3);
  Pi[56]=P(ix3,iy2,iz0);  Pi[57]=P(ix3,iy2,iz1);  Pi[58]=P(ix3,iy2,iz2);  Pi[59]=P(ix3,iy2,iz3);
  Pi[60]=P(ix3,iy3,iz0);  Pi[61]=P(ix3,iy3,iz1);  Pi[62]=P(ix3,iy3,iz2);  Pi[63]=P(ix3,iy3,iz3);

  cP[ 0] = c[0]*Pi[ 0] + c[1]*Pi[ 1] + c[2]*Pi[ 2] + c[3]*Pi[ 3];
  cP[ 1] = c[0]*Pi[ 4] + c[1]*Pi[ 5] + c[2]*Pi[ 6] + c[3]*Pi[ 7];
  cP[ 2] = c[0]*Pi[ 8] + c[1]*Pi[ 9] + c[2]*Pi[10] + c[3]*Pi[11];
  cP[ 3] = c[0]*Pi[12] + c[1]*Pi[13] + c[2]*Pi[14] + c[3]*Pi[15];
  cP[ 4] = c[0]*Pi[16] + c[1]*Pi[17] + c[2]*Pi[18] + c[3]*Pi[19];
  cP[ 5] = c[0]*Pi[20] + c[1]*Pi[21] + c[2]*Pi[22] + c[3]*Pi[23];
  cP[ 6] = c[0]*Pi[24] + c[1]*Pi[25] + c[2]*Pi[26] + c[3]*Pi[27];
  cP[ 7] = c[0]*Pi[28] + c[1]*Pi[29] + c[2]*Pi[30] + c[3]*Pi[31];
  cP[ 8] = c[0]*Pi[32] + c[1]*Pi[33] + c[2]*Pi[34] + c[3]*Pi[35];
  cP[ 9] = c[0]*Pi[36] + c[1]*Pi[37] + c[2]*Pi[38] + c[3]*Pi[39];
  cP[10] = c[0]*Pi[40] + c[1]*Pi[41] + c[2]*Pi[42] + c[3]*Pi[43];
  cP[11] = c[0]*Pi[44] + c[1]*Pi[45] + c[2]*Pi[46] + c[3]*Pi[47];
  cP[12] = c[0]*Pi[48] + c[1]*Pi[49] + c[2]*Pi[50] + c[3]*Pi[51];
  cP[13] = c[0]*Pi[52] + c[1]*Pi[53] + c[2]*Pi[54] + c[3]*Pi[55];
  cP[14] = c[0]*Pi[56] + c[1]*Pi[57] + c[2]*Pi[58] + c[3]*Pi[59];
  cP[15] = c[0]*Pi[60] + c[1]*Pi[61] + c[2]*Pi[62] + c[3]*Pi[63];

  bcP[0] = b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3];
  bcP[1] = b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7];
  bcP[2] = b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11];
  bcP[3] = b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15];

  val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  grad[0] = da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3];
  grad[1] = (a[0]*(db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]) +
	     a[1]*(db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]) +
	     a[2]*(db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]) +
	     a[3]*(db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]));
  grad[2] = 
    a[0]*(b[0]*(dc[0]*Pi[ 0] + dc[1]*Pi[ 1] + dc[2]*Pi[ 2] + dc[3]*Pi[ 3])+
	  b[1]*(dc[0]*Pi[ 4] + dc[1]*Pi[ 5] + dc[2]*Pi[ 6] + dc[3]*Pi[ 7])+
	  b[2]*(dc[0]*Pi[ 8] + dc[1]*Pi[ 9] + dc[2]*Pi[10] + dc[3]*Pi[11])+
	  b[3]*(dc[0]*Pi[12] + dc[1]*Pi[13] + dc[2]*Pi[14] + dc[3]*Pi[15]))+
    a[1]*(b[0]*(dc[0]*Pi[16] + dc[1]*Pi[17] + dc[2]*Pi[18] + dc[3]*Pi[19])+ 
	  b[1]*(dc[0]*Pi[20] + dc[1]*Pi[21] + dc[2]*Pi[22] + dc[3]*Pi[23])+
	  b[2]*(dc[0]*Pi[24] + dc[1]*Pi[25] + dc[2]*Pi[26] + dc[3]*Pi[27])+
	  b[3]*(dc[0]*Pi[28] + dc[1]*Pi[29] + dc[2]*Pi[30] + dc[3]*Pi[31]))+
    a[2]*(b[0]*(dc[0]*Pi[32] + dc[1]*Pi[33] + dc[2]*Pi[34] + dc[3]*Pi[35])+
	  b[1]*(dc[0]*Pi[36] + dc[1]*Pi[37] + dc[2]*Pi[38] + dc[3]*Pi[39])+
	  b[2]*(dc[0]*Pi[40] + dc[1]*Pi[41] + dc[2]*Pi[42] + dc[3]*Pi[43])+
	  b[3]*(dc[0]*Pi[44] + dc[1]*Pi[45] + dc[2]*Pi[46] + dc[3]*Pi[47]))+
    a[3]*(b[0]*(dc[0]*Pi[48] + dc[1]*Pi[49] + dc[2]*Pi[50] + dc[3]*Pi[51])+
	  b[1]*(dc[0]*Pi[52] + dc[1]*Pi[53] + dc[2]*Pi[54] + dc[3]*Pi[55])+
	  b[2]*(dc[0]*Pi[56] + dc[1]*Pi[57] + dc[2]*Pi[58] + dc[3]*Pi[59])+
	  b[3]*(dc[0]*Pi[60] + dc[1]*Pi[61] + dc[2]*Pi[62] + dc[3]*Pi[63]));
}


template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Evaluate (TinyVector<double,3> r,
							      T &val,
							      TinyVector<T,3> &grad,
							      T &laplacian) const
{
  TinyVector<double,4> a, da, d2a, b, db, d2b, c, dc, d2c;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a, da, d2a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b, db, d2b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c, dc, d2c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  T Pi[64], cP[16], bcP[4];
  Pi[ 0]=P(ix0,iy0,iz0);  Pi[ 1]=P(ix0,iy0,iz1);  Pi[ 2]=P(ix0,iy0,iz2);  Pi[ 3]=P(ix0,iy0,iz3);
  Pi[ 4]=P(ix0,iy1,iz0);  Pi[ 5]=P(ix0,iy1,iz1);  Pi[ 6]=P(ix0,iy1,iz2);  Pi[ 7]=P(ix0,iy1,iz3);
  Pi[ 8]=P(ix0,iy2,iz0);  Pi[ 9]=P(ix0,iy2,iz1);  Pi[10]=P(ix0,iy2,iz2);  Pi[11]=P(ix0,iy2,iz3);
  Pi[12]=P(ix0,iy3,iz0);  Pi[13]=P(ix0,iy3,iz1);  Pi[14]=P(ix0,iy3,iz2);  Pi[15]=P(ix0,iy3,iz3);
  Pi[16]=P(ix1,iy0,iz0);  Pi[17]=P(ix1,iy0,iz1);  Pi[18]=P(ix1,iy0,iz2);  Pi[19]=P(ix1,iy0,iz3);
  Pi[20]=P(ix1,iy1,iz0);  Pi[21]=P(ix1,iy1,iz1);  Pi[22]=P(ix1,iy1,iz2);  Pi[23]=P(ix1,iy1,iz3);
  Pi[24]=P(ix1,iy2,iz0);  Pi[25]=P(ix1,iy2,iz1);  Pi[26]=P(ix1,iy2,iz2);  Pi[27]=P(ix1,iy2,iz3);
  Pi[28]=P(ix1,iy3,iz0);  Pi[29]=P(ix1,iy3,iz1);  Pi[30]=P(ix1,iy3,iz2);  Pi[31]=P(ix1,iy3,iz3);
  Pi[32]=P(ix2,iy0,iz0);  Pi[33]=P(ix2,iy0,iz1);  Pi[34]=P(ix2,iy0,iz2);  Pi[35]=P(ix2,iy0,iz3);
  Pi[36]=P(ix2,iy1,iz0);  Pi[37]=P(ix2,iy1,iz1);  Pi[38]=P(ix2,iy1,iz2);  Pi[39]=P(ix2,iy1,iz3);
  Pi[40]=P(ix2,iy2,iz0);  Pi[41]=P(ix2,iy2,iz1);  Pi[42]=P(ix2,iy2,iz2);  Pi[43]=P(ix2,iy2,iz3);
  Pi[44]=P(ix2,iy3,iz0);  Pi[45]=P(ix2,iy3,iz1);  Pi[46]=P(ix2,iy3,iz2);  Pi[47]=P(ix2,iy3,iz3);
  Pi[48]=P(ix3,iy0,iz0);  Pi[49]=P(ix3,iy0,iz1);  Pi[50]=P(ix3,iy0,iz2);  Pi[51]=P(ix3,iy0,iz3);
  Pi[52]=P(ix3,iy1,iz0);  Pi[53]=P(ix3,iy1,iz1);  Pi[54]=P(ix3,iy1,iz2);  Pi[55]=P(ix3,iy1,iz3);
  Pi[56]=P(ix3,iy2,iz0);  Pi[57]=P(ix3,iy2,iz1);  Pi[58]=P(ix3,iy2,iz2);  Pi[59]=P(ix3,iy2,iz3);
  Pi[60]=P(ix3,iy3,iz0);  Pi[61]=P(ix3,iy3,iz1);  Pi[62]=P(ix3,iy3,iz2);  Pi[63]=P(ix3,iy3,iz3);

  cP[ 0] = c[0]*Pi[ 0] + c[1]*Pi[ 1] + c[2]*Pi[ 2] + c[3]*Pi[ 3];
  cP[ 1] = c[0]*Pi[ 4] + c[1]*Pi[ 5] + c[2]*Pi[ 6] + c[3]*Pi[ 7];
  cP[ 2] = c[0]*Pi[ 8] + c[1]*Pi[ 9] + c[2]*Pi[10] + c[3]*Pi[11];
  cP[ 3] = c[0]*Pi[12] + c[1]*Pi[13] + c[2]*Pi[14] + c[3]*Pi[15];
  cP[ 4] = c[0]*Pi[16] + c[1]*Pi[17] + c[2]*Pi[18] + c[3]*Pi[19];
  cP[ 5] = c[0]*Pi[20] + c[1]*Pi[21] + c[2]*Pi[22] + c[3]*Pi[23];
  cP[ 6] = c[0]*Pi[24] + c[1]*Pi[25] + c[2]*Pi[26] + c[3]*Pi[27];
  cP[ 7] = c[0]*Pi[28] + c[1]*Pi[29] + c[2]*Pi[30] + c[3]*Pi[31];
  cP[ 8] = c[0]*Pi[32] + c[1]*Pi[33] + c[2]*Pi[34] + c[3]*Pi[35];
  cP[ 9] = c[0]*Pi[36] + c[1]*Pi[37] + c[2]*Pi[38] + c[3]*Pi[39];
  cP[10] = c[0]*Pi[40] + c[1]*Pi[41] + c[2]*Pi[42] + c[3]*Pi[43];
  cP[11] = c[0]*Pi[44] + c[1]*Pi[45] + c[2]*Pi[46] + c[3]*Pi[47];
  cP[12] = c[0]*Pi[48] + c[1]*Pi[49] + c[2]*Pi[50] + c[3]*Pi[51];
  cP[13] = c[0]*Pi[52] + c[1]*Pi[53] + c[2]*Pi[54] + c[3]*Pi[55];
  cP[14] = c[0]*Pi[56] + c[1]*Pi[57] + c[2]*Pi[58] + c[3]*Pi[59];
  cP[15] = c[0]*Pi[60] + c[1]*Pi[61] + c[2]*Pi[62] + c[3]*Pi[63];

  bcP[0] = b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3];
  bcP[1] = b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7];
  bcP[2] = b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11];
  bcP[3] = b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15];

  val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  grad[0] = da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3];
  grad[1] = (a[0]*(db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]) +
	     a[1]*(db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]) +
	     a[2]*(db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]) +
	     a[3]*(db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]));
  grad[2] = 
    a[0]*(b[0]*(dc[0]*Pi[ 0] + dc[1]*Pi[ 1] + dc[2]*Pi[ 2] + dc[3]*Pi[ 3])+
	  b[1]*(dc[0]*Pi[ 4] + dc[1]*Pi[ 5] + dc[2]*Pi[ 6] + dc[3]*Pi[ 7])+
	  b[2]*(dc[0]*Pi[ 8] + dc[1]*Pi[ 9] + dc[2]*Pi[10] + dc[3]*Pi[11])+
	  b[3]*(dc[0]*Pi[12] + dc[1]*Pi[13] + dc[2]*Pi[14] + dc[3]*Pi[15]))+
    a[1]*(b[0]*(dc[0]*Pi[16] + dc[1]*Pi[17] + dc[2]*Pi[18] + dc[3]*Pi[19])+ 
	  b[1]*(dc[0]*Pi[20] + dc[1]*Pi[21] + dc[2]*Pi[22] + dc[3]*Pi[23])+
	  b[2]*(dc[0]*Pi[24] + dc[1]*Pi[25] + dc[2]*Pi[26] + dc[3]*Pi[27])+
	  b[3]*(dc[0]*Pi[28] + dc[1]*Pi[29] + dc[2]*Pi[30] + dc[3]*Pi[31]))+
    a[2]*(b[0]*(dc[0]*Pi[32] + dc[1]*Pi[33] + dc[2]*Pi[34] + dc[3]*Pi[35])+
	  b[1]*(dc[0]*Pi[36] + dc[1]*Pi[37] + dc[2]*Pi[38] + dc[3]*Pi[39])+
	  b[2]*(dc[0]*Pi[40] + dc[1]*Pi[41] + dc[2]*Pi[42] + dc[3]*Pi[43])+
	  b[3]*(dc[0]*Pi[44] + dc[1]*Pi[45] + dc[2]*Pi[46] + dc[3]*Pi[47]))+
    a[3]*(b[0]*(dc[0]*Pi[48] + dc[1]*Pi[49] + dc[2]*Pi[50] + dc[3]*Pi[51])+
	  b[1]*(dc[0]*Pi[52] + dc[1]*Pi[53] + dc[2]*Pi[54] + dc[3]*Pi[55])+
	  b[2]*(dc[0]*Pi[56] + dc[1]*Pi[57] + dc[2]*Pi[58] + dc[3]*Pi[59])+
	  b[3]*(dc[0]*Pi[60] + dc[1]*Pi[61] + dc[2]*Pi[62] + dc[3]*Pi[63]));

  laplacian = d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3];
  laplacian += (a[0]*(d2b[0]*cP[ 0] + d2b[1]*cP[ 1] + d2b[2]*cP[ 2] + d2b[3]*cP[ 3]) +
		a[1]*(d2b[0]*cP[ 4] + d2b[1]*cP[ 5] + d2b[2]*cP[ 6] + d2b[3]*cP[ 7]) +
		a[2]*(d2b[0]*cP[ 8] + d2b[1]*cP[ 9] + d2b[2]*cP[10] + d2b[3]*cP[11]) +
		a[3]*(d2b[0]*cP[12] + d2b[1]*cP[13] + d2b[2]*cP[14] + d2b[3]*cP[15]));
  laplacian +=     
    a[0]*(b[0]*(d2c[0]*Pi[ 0] + d2c[1]*Pi[ 1] + d2c[2]*Pi[ 2] + d2c[3]*Pi[ 3])+
	  b[1]*(d2c[0]*Pi[ 4] + d2c[1]*Pi[ 5] + d2c[2]*Pi[ 6] + d2c[3]*Pi[ 7])+
	  b[2]*(d2c[0]*Pi[ 8] + d2c[1]*Pi[ 9] + d2c[2]*Pi[10] + d2c[3]*Pi[11])+
	  b[3]*(d2c[0]*Pi[12] + d2c[1]*Pi[13] + d2c[2]*Pi[14] + d2c[3]*Pi[15]))+
    a[1]*(b[0]*(d2c[0]*Pi[16] + d2c[1]*Pi[17] + d2c[2]*Pi[18] + d2c[3]*Pi[19])+ 
	  b[1]*(d2c[0]*Pi[20] + d2c[1]*Pi[21] + d2c[2]*Pi[22] + d2c[3]*Pi[23])+
	  b[2]*(d2c[0]*Pi[24] + d2c[1]*Pi[25] + d2c[2]*Pi[26] + d2c[3]*Pi[27])+
	  b[3]*(d2c[0]*Pi[28] + d2c[1]*Pi[29] + d2c[2]*Pi[30] + d2c[3]*Pi[31]))+
    a[2]*(b[0]*(d2c[0]*Pi[32] + d2c[1]*Pi[33] + d2c[2]*Pi[34] + d2c[3]*Pi[35])+
	  b[1]*(d2c[0]*Pi[36] + d2c[1]*Pi[37] + d2c[2]*Pi[38] + d2c[3]*Pi[39])+
	  b[2]*(d2c[0]*Pi[40] + d2c[1]*Pi[41] + d2c[2]*Pi[42] + d2c[3]*Pi[43])+
	  b[3]*(d2c[0]*Pi[44] + d2c[1]*Pi[45] + d2c[2]*Pi[46] + d2c[3]*Pi[47]))+
    a[3]*(b[0]*(d2c[0]*Pi[48] + d2c[1]*Pi[49] + d2c[2]*Pi[50] + d2c[3]*Pi[51])+
	  b[1]*(d2c[0]*Pi[52] + d2c[1]*Pi[53] + d2c[2]*Pi[54] + d2c[3]*Pi[55])+
	  b[2]*(d2c[0]*Pi[56] + d2c[1]*Pi[57] + d2c[2]*Pi[58] + d2c[3]*Pi[59])+
	  b[3]*(d2c[0]*Pi[60] + d2c[1]*Pi[61] + d2c[2]*Pi[62] + d2c[3]*Pi[63]));

}




template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Evaluate 
(TinyVector<double,3> r, T &val, TinyVector<T,3> &grad, TinyMatrix<T,3,3> &secDerivs) const
{
  TinyVector<double,4> a, da, d2a, b, db, d2b, c, dc, d2c;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a, da, d2a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b, db, d2b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c, dc, d2c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  // Save some operations by factorizing computation.
  TinyMatrix<T,4,4> cP, dcP;
  cP(0,0) = c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3);
  cP(0,1) = c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3);
  cP(0,2) = c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3);
  cP(0,3) = c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3);
  cP(1,0) = c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3);
  cP(1,1) = c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3); 
  cP(1,2) = c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3);
  cP(1,3) = c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3);
  cP(2,0) = c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3);
  cP(2,1) = c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3);
  cP(2,2) = c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3);
  cP(2,3) = c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3);
  cP(3,0) = c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3);
  cP(3,1) = c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3);
  cP(3,2) = c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3);
  cP(3,3) = c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3);

  dcP(0,0) = dc[0]*P(ix0,iy0,iz0)+dc[1]*P(ix0,iy0,iz1)+dc[2]*P(ix0,iy0,iz2)+dc[3]*P(ix0,iy0,iz3);
  dcP(0,1) = dc[0]*P(ix0,iy1,iz0)+dc[1]*P(ix0,iy1,iz1)+dc[2]*P(ix0,iy1,iz2)+dc[3]*P(ix0,iy1,iz3);
  dcP(0,2) = dc[0]*P(ix0,iy2,iz0)+dc[1]*P(ix0,iy2,iz1)+dc[2]*P(ix0,iy2,iz2)+dc[3]*P(ix0,iy2,iz3);
  dcP(0,3) = dc[0]*P(ix0,iy3,iz0)+dc[1]*P(ix0,iy3,iz1)+dc[2]*P(ix0,iy3,iz2)+dc[3]*P(ix0,iy3,iz3);
  dcP(1,0) = dc[0]*P(ix1,iy0,iz0)+dc[1]*P(ix1,iy0,iz1)+dc[2]*P(ix1,iy0,iz2)+dc[3]*P(ix1,iy0,iz3);
  dcP(1,1) = dc[0]*P(ix1,iy1,iz0)+dc[1]*P(ix1,iy1,iz1)+dc[2]*P(ix1,iy1,iz2)+dc[3]*P(ix1,iy1,iz3); 
  dcP(1,2) = dc[0]*P(ix1,iy2,iz0)+dc[1]*P(ix1,iy2,iz1)+dc[2]*P(ix1,iy2,iz2)+dc[3]*P(ix1,iy2,iz3);
  dcP(1,3) = dc[0]*P(ix1,iy3,iz0)+dc[1]*P(ix1,iy3,iz1)+dc[2]*P(ix1,iy3,iz2)+dc[3]*P(ix1,iy3,iz3);
  dcP(2,0) = dc[0]*P(ix2,iy0,iz0)+dc[1]*P(ix2,iy0,iz1)+dc[2]*P(ix2,iy0,iz2)+dc[3]*P(ix2,iy0,iz3);
  dcP(2,1) = dc[0]*P(ix2,iy1,iz0)+dc[1]*P(ix2,iy1,iz1)+dc[2]*P(ix2,iy1,iz2)+dc[3]*P(ix2,iy1,iz3);
  dcP(2,2) = dc[0]*P(ix2,iy2,iz0)+dc[1]*P(ix2,iy2,iz1)+dc[2]*P(ix2,iy2,iz2)+dc[3]*P(ix2,iy2,iz3);
  dcP(2,3) = dc[0]*P(ix2,iy3,iz0)+dc[1]*P(ix2,iy3,iz1)+dc[2]*P(ix2,iy3,iz2)+dc[3]*P(ix2,iy3,iz3);
  dcP(3,0) = dc[0]*P(ix3,iy0,iz0)+dc[1]*P(ix3,iy0,iz1)+dc[2]*P(ix3,iy0,iz2)+dc[3]*P(ix3,iy0,iz3);
  dcP(3,1) = dc[0]*P(ix3,iy1,iz0)+dc[1]*P(ix3,iy1,iz1)+dc[2]*P(ix3,iy1,iz2)+dc[3]*P(ix3,iy1,iz3);
  dcP(3,2) = dc[0]*P(ix3,iy2,iz0)+dc[1]*P(ix3,iy2,iz1)+dc[2]*P(ix3,iy2,iz2)+dc[3]*P(ix3,iy2,iz3);
  dcP(3,3) = dc[0]*P(ix3,iy3,iz0)+dc[1]*P(ix3,iy3,iz1)+dc[2]*P(ix3,iy3,iz2)+dc[3]*P(ix3,iy3,iz3);

  TinyVector<T,4> bcP, bdcP;
  bcP[0] = cP(0,0)*b[0]+cP(0,1)*b[1]+cP(0,2)*b[2]+cP(0,3)*b[3];
  bcP[1] = cP(1,0)*b[0]+cP(1,1)*b[1]+cP(1,2)*b[2]+cP(1,3)*b[3];
  bcP[2] = cP(2,0)*b[0]+cP(2,1)*b[1]+cP(2,2)*b[2]+cP(2,3)*b[3];
  bcP[3] = cP(3,0)*b[0]+cP(3,1)*b[1]+cP(3,2)*b[2]+cP(3,3)*b[3];
 
  bdcP[0] = dcP(0,0)*b[0]+dcP(0,1)*b[1]+dcP(0,2)*b[2]+dcP(0,3)*b[3];
  bdcP[1] = dcP(1,0)*b[0]+dcP(1,1)*b[1]+dcP(1,2)*b[2]+dcP(1,3)*b[3];
  bdcP[2] = dcP(2,0)*b[0]+dcP(2,1)*b[1]+dcP(2,2)*b[2]+dcP(2,3)*b[3];
  bdcP[3] = dcP(3,0)*b[0]+dcP(3,1)*b[1]+dcP(3,2)*b[2]+dcP(3,3)*b[3];



  // Compute value
  val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];

  // Compute gradient
  grad[0] = 
    (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  
  grad[1] =
    (a[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
     a[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
     a[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
     a[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));

  grad[2] = 
    (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);


  // Compute laplacian
  secDerivs(0,0) = 
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3]);
  secDerivs(0,1) = secDerivs(1,0) = 
    (da[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
     da[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
     da[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
     da[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));
  secDerivs(0,2) = secDerivs(2,0) = 
    (da[0]*bdcP[0] + da[1]*bdcP[1] + da[2]*bdcP[2] + da[3]*bdcP[3]);
        
  secDerivs(1,1) = 
    (a[0]*(cP(0,0)*d2b[0]+cP(0,1)*d2b[1]+cP(0,2)*d2b[2]+cP(0,3)*d2b[3]) +
     a[1]*(cP(1,0)*d2b[0]+cP(1,1)*d2b[1]+cP(1,2)*d2b[2]+cP(1,3)*d2b[3]) +
     a[2]*(cP(2,0)*d2b[0]+cP(2,1)*d2b[1]+cP(2,2)*d2b[2]+cP(2,3)*d2b[3]) +
     a[3]*(cP(3,0)*d2b[0]+cP(3,1)*d2b[1]+cP(3,2)*d2b[2]+cP(3,3)*d2b[3]));
  secDerivs(1,2) = secDerivs(2,1) = 
    (a[0]*(db[0]*dcP(0,0) + db[1]*dcP(0,1) + db[2]*dcP(0,2) + db[3]*dcP(0,3))+
     a[1]*(db[0]*dcP(1,0) + db[1]*dcP(1,1) + db[2]*dcP(1,2) + db[3]*dcP(1,3))+
     a[2]*(db[0]*dcP(2,0) + db[1]*dcP(2,1) + db[2]*dcP(2,2) + db[3]*dcP(2,3))+
     a[3]*(db[0]*dcP(3,0) + db[1]*dcP(3,1) + db[2]*dcP(3,2) + db[3]*dcP(3,3)));
  secDerivs(2,2) = 
    (a[0]*(b[0]*(d2c[0]*P(ix0,iy0,iz0)+d2c[1]*P(ix0,iy0,iz1)+d2c[2]*P(ix0,iy0,iz2)+d2c[3]*P(ix0,iy0,iz3))+
	   b[1]*(d2c[0]*P(ix0,iy1,iz0)+d2c[1]*P(ix0,iy1,iz1)+d2c[2]*P(ix0,iy1,iz2)+d2c[3]*P(ix0,iy1,iz3))+
	   b[2]*(d2c[0]*P(ix0,iy2,iz0)+d2c[1]*P(ix0,iy2,iz1)+d2c[2]*P(ix0,iy2,iz2)+d2c[3]*P(ix0,iy2,iz3))+
	   b[3]*(d2c[0]*P(ix0,iy3,iz0)+d2c[1]*P(ix0,iy3,iz1)+d2c[2]*P(ix0,iy3,iz2)+d2c[3]*P(ix0,iy3,iz3)))+
     a[1]*(b[0]*(d2c[0]*P(ix1,iy0,iz0)+d2c[1]*P(ix1,iy0,iz1)+d2c[2]*P(ix1,iy0,iz2)+d2c[3]*P(ix1,iy0,iz3))+
	   b[1]*(d2c[0]*P(ix1,iy1,iz0)+d2c[1]*P(ix1,iy1,iz1)+d2c[2]*P(ix1,iy1,iz2)+d2c[3]*P(ix1,iy1,iz3))+
	   b[2]*(d2c[0]*P(ix1,iy2,iz0)+d2c[1]*P(ix1,iy2,iz1)+d2c[2]*P(ix1,iy2,iz2)+d2c[3]*P(ix1,iy2,iz3))+
	   b[3]*(d2c[0]*P(ix1,iy3,iz0)+d2c[1]*P(ix1,iy3,iz1)+d2c[2]*P(ix1,iy3,iz2)+d2c[3]*P(ix1,iy3,iz3)))+
     a[2]*(b[0]*(d2c[0]*P(ix2,iy0,iz0)+d2c[1]*P(ix2,iy0,iz1)+d2c[2]*P(ix2,iy0,iz2)+d2c[3]*P(ix2,iy0,iz3))+
	   b[1]*(d2c[0]*P(ix2,iy1,iz0)+d2c[1]*P(ix2,iy1,iz1)+d2c[2]*P(ix2,iy1,iz2)+d2c[3]*P(ix2,iy1,iz3))+
	   b[2]*(d2c[0]*P(ix2,iy2,iz0)+d2c[1]*P(ix2,iy2,iz1)+d2c[2]*P(ix2,iy2,iz2)+d2c[3]*P(ix2,iy2,iz3))+
	   b[3]*(d2c[0]*P(ix2,iy3,iz0)+d2c[1]*P(ix2,iy3,iz1)+d2c[2]*P(ix2,iy3,iz2)+d2c[3]*P(ix2,iy3,iz3)))+
     a[3]*(b[0]*(d2c[0]*P(ix3,iy0,iz0)+d2c[1]*P(ix3,iy0,iz1)+d2c[2]*P(ix3,iy0,iz2)+d2c[3]*P(ix3,iy0,iz3))+
	   b[1]*(d2c[0]*P(ix3,iy1,iz0)+d2c[1]*P(ix3,iy1,iz1)+d2c[2]*P(ix3,iy1,iz2)+d2c[3]*P(ix3,iy1,iz3))+
	   b[2]*(d2c[0]*P(ix3,iy2,iz0)+d2c[1]*P(ix3,iy2,iz1)+d2c[2]*P(ix3,iy2,iz2)+d2c[3]*P(ix3,iy2,iz3))+
	   b[3]*(d2c[0]*P(ix3,iy3,iz0)+d2c[1]*P(ix3,iy3,iz1)+d2c[2]*P(ix3,iy3,iz2)+d2c[3]*P(ix3,iy3,iz3))));

}


//////////////////////////////////////////////////////////////////////
//           Optimized, SSE version for floating point              //
//////////////////////////////////////////////////////////////////////
#ifdef __SSE3__
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

template<typename XGridType, typename YGridType, typename ZGridType>
class TricubicNUBspline<float,XGridType,YGridType,ZGridType>
{
private:
  NUBsplineBasis<XGridType> XBasis;
  NUBsplineBasis<YGridType> YBasis;
  NUBsplineBasis<ZGridType> ZBasis;
  int Nx, Ny, Nz;
  // The starting and ending values for the uniform grids
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  // The box dimensions and their inverses
  double Lx, LxInv, Ly, LyInv, Lz, LzInv;
  // The control points
  Array<float,3> P;

  TinyVector<double,3> Periodic;  
public:
  inline void Init (XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<float,3> &data, 
		    BCType xbct=PERIODIC, BCType ybct=PERIODIC, BCType zbct=PERIODIC);
  inline float operator() (TinyVector<double,3> r) const;
  inline void Evaluate (TinyVector<double,3> r, float &val, 
			TinyVector<float,3> &grad) const;
  inline void Evaluate (TinyVector<double,3> r, float &val,
			TinyVector<float,3> &grad, float &laplacian) const;
  inline void Evaluate (TinyVector<double,3> r, float & val,
			TinyVector<float,3> &grad, 
			TinyMatrix<float,3,3> &secDerivs) const;
};

template<typename XGridType, typename YGridType, typename ZGridType>
void 
TricubicNUBspline<float,XGridType,YGridType,ZGridType>::Init
(XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<float,3> &data, 
 BCType xbc, BCType ybc, BCType zbc)
{
  // Set up 1D basis functions
  XBasis.Init (xgrid, xbc==PERIODIC);
  YBasis.Init (ygrid, ybc==PERIODIC);
  ZBasis.Init (zgrid, zbc==PERIODIC);
  Periodic[0] = xbc==PERIODIC ? 1.0 : 0.0;
  Periodic[1] = ybc==PERIODIC ? 1.0 : 0.0;
  Periodic[2] = zbc==PERIODIC ? 1.0 : 0.0;

  Nx = data.extent(0);
  Ny = data.extent(1);
  Nz = data.extent(2);

  int Mx, My, Mz;
  Mx = (xbc == PERIODIC) ? Nx+3 : Nx+2;
  My = (ybc == PERIODIC) ? Ny+3 : Ny+2;
  Mz = (zbc == PERIODIC) ? Nz+3 : Nz+2;
  assert (xgrid->NumPoints == Mx-2);
  assert (ygrid->NumPoints == My-2);
  assert (zgrid->NumPoints == Mz-2);

  P.resize(Mx, My, Mz);

  
  // Now solve interpolating equations
  TinyVector<float,4> lBC, rBC, dummy1, dummy2;
  ////////////////////
  // Do X direction //
  ////////////////////
  if (xbc == PERIODIC) {
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	SolvePeriodicInterp1D(XBasis, data(Range::all(), iy, iz), P(Range::all(),iy+1, iz+1));
  }
  else {
    if (xbc == FLAT) {
      XBasis (0   , dummy1, lBC);
      XBasis (Nx-1, dummy1, rBC);
    }
    else if (xbc == NATURAL) {
      XBasis (0   , dummy1, dummy2, lBC);
      XBasis (Nx-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	SolveDerivInterp1D(XBasis, data(Range::all(), iy, iz),  P(Range::all(), iy+1, iz+1), lBC, rBC);
  }
  ////////////////////
  // Do Y direction //
  ////////////////////
  if (ybc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++) 
	SolvePeriodicInterp1D(YBasis, P(ix, Range(1,Ny), iz), P(ix, Range::all(), iz));
  }
  else {
    if (ybc == FLAT) {
      YBasis (0,    dummy1, lBC);
      YBasis (Ny-1, dummy1, rBC);
    }
    else if (ybc == NATURAL) {
      YBasis (0,    dummy1, dummy2, lBC);
      YBasis (Ny-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++)
	SolveDerivInterp1D(YBasis, P(ix,Range(1,Ny), iz),  P(ix, Range::all(), iz), lBC, rBC);
  }
  ////////////////////
  // Do Z direction //
  ////////////////////
  if (zbc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++) 
	SolvePeriodicInterp1D(ZBasis, P(ix, iy, Range(1,Nz)), P(ix, iy, Range::all()));
  }
  else {
    if (zbc == FLAT) {
      ZBasis (0,    dummy1, lBC);
      ZBasis (Nz-1, dummy1, rBC);
    }
    else if (zbc == NATURAL) {
      ZBasis (0,    dummy1, dummy2, lBC);
      ZBasis (Nz-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++)
	SolveDerivInterp1D(ZBasis, P(ix, iy, Range(1,Nz)),  P(ix, iy, Range::all()), lBC, rBC);
  }
}


template<typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<float,XGridType,YGridType,ZGridType>::Evaluate (TinyVector<double,3> r,
								  float &val,
								  TinyVector<float,3> &grad,
								  float &laplacian) const
{
  TinyVector<float,4> a, da, d2a, b, db, d2b, c, dc, d2c;
  __m128 av, dav, d2av, bv, dbv, d2bv, cv, dcv, d2cv;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a, da, d2a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b, db, d2b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c, dc, d2c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;
  av = _mm_loadu_ps (&(a[0]));  dav = _mm_loadu_ps (&(da[0])); d2av = _mm_loadu_ps (&(d2a[0]));
  bv = _mm_loadu_ps (&(b[0]));  dbv = _mm_loadu_ps (&(db[0])); d2bv = _mm_loadu_ps (&(d2b[0]));
  cv = _mm_loadu_ps (&(c[0]));  dcv = _mm_loadu_ps (&(dc[0])); d2cv = _mm_loadu_ps (&(d2c[0]));

  __m128 Pi[16], cP[4], bcP;
  Pi[0]  = _mm_loadu_ps (&P(ix0, iy0, iz0));
  Pi[1]  = _mm_loadu_ps (&P(ix0, iy1, iz0));
  Pi[2]  = _mm_loadu_ps (&P(ix0, iy2, iz0));
  Pi[3]  = _mm_loadu_ps (&P(ix0, iy3, iz0));
  Pi[4]  = _mm_loadu_ps (&P(ix1, iy0, iz0));
  Pi[5]  = _mm_loadu_ps (&P(ix1, iy1, iz0));
  Pi[6]  = _mm_loadu_ps (&P(ix1, iy2, iz0));
  Pi[7]  = _mm_loadu_ps (&P(ix1, iy3, iz0));
  Pi[8]  = _mm_loadu_ps (&P(ix2, iy0, iz0));
  Pi[9]  = _mm_loadu_ps (&P(ix2, iy1, iz0));
  Pi[10] = _mm_loadu_ps (&P(ix2, iy2, iz0));
  Pi[11] = _mm_loadu_ps (&P(ix2, iy3, iz0));
  Pi[12] = _mm_loadu_ps (&P(ix3, iy0, iz0));
  Pi[13] = _mm_loadu_ps (&P(ix3, iy1, iz0));
  Pi[14] = _mm_loadu_ps (&P(ix3, iy2, iz0));
  Pi[15] = _mm_loadu_ps (&P(ix3, iy3, iz0));


}

template<typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<float,XGridType,YGridType,ZGridType>::Evaluate 
(TinyVector<double,3> r, float &val, TinyVector<float,3> &grad, TinyMatrix<float,3,3> &secDerivs) const
{
  TinyVector<float,4>  af, daf, d2af, bf, dbf, d2bf, cf, dcf, d2cf;
  __m128 av, dav, d2av, bv, dbv, d2bv, cv, dcv, d2cv, reg0, reg1, reg2, reg3;
  // Evaluate 1D basis functions
//   int ix0 = XBasis(r[0], a, da, d2a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
//   int iy0 = YBasis(r[1], b, db, d2b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
//   int iz0 = ZBasis(r[2], c, dc, d2c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;
//   af=a; daf=da; d2af = d2a;
//   bf=b; dbf=db; d2bf = d2b;
//   cf=c; dcf=dc; d2cf = d2c;
   int ix0 = XBasis(r[0], af, daf, d2af); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
   int iy0 = YBasis(r[1], bf, dbf, d2bf); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
   int iz0 = ZBasis(r[2], cf, dcf, d2cf); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  av = _mm_loadu_ps (&(af[0]));  dav = _mm_loadu_ps (&(daf[0])); d2av = _mm_loadu_ps (&(d2af[0]));
  bv = _mm_loadu_ps (&(bf[0]));  dbv = _mm_loadu_ps (&(dbf[0])); d2bv = _mm_loadu_ps (&(d2bf[0]));
  cv = _mm_loadu_ps (&(cf[0]));  dcv = _mm_loadu_ps (&(dcf[0])); d2cv = _mm_loadu_ps (&(d2cf[0]));

  __m128 tmp[16], cP[4], dcP[4], d2cP[4], bcP, bdcP;
  reg0 = _mm_loadu_ps (&(P(ix0, iy0, iz0)));    
  reg1 = _mm_loadu_ps (&(P(ix0, iy1, iz0)));  
  reg2 = _mm_loadu_ps (&(P(ix0, iy2, iz0)));  
  reg3 = _mm_loadu_ps (&(P(ix0, iy3, iz0)));  
  // cP
  tmp[0]   = _mm_mul_ps (cv, reg0);
  tmp[1]   = _mm_mul_ps (cv, reg1);
  tmp[2]   = _mm_mul_ps (cv, reg2);
  tmp[3]   = _mm_mul_ps (cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  cP[0]    = _mm_hadd_ps (tmp[0], tmp[1]);
  // dcP
  tmp[0]   = _mm_mul_ps (dcv, reg0);
  tmp[1]   = _mm_mul_ps (dcv, reg1);
  tmp[2]   = _mm_mul_ps (dcv, reg2);
  tmp[3]   = _mm_mul_ps (dcv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  dcP[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  // d2cP
  tmp[0]   = _mm_mul_ps (d2cv, reg0);
  tmp[1]   = _mm_mul_ps (d2cv, reg1);
  tmp[2]   = _mm_mul_ps (d2cv, reg2);
  tmp[3]   = _mm_mul_ps (d2cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  d2cP[0]  = _mm_hadd_ps (tmp[0], tmp[1]);

  reg0 = _mm_loadu_ps (&(P(ix1, iy0, iz0)));    
  reg1 = _mm_loadu_ps (&(P(ix1, iy1, iz0)));  
  reg2 = _mm_loadu_ps (&(P(ix1, iy2, iz0)));  
  reg3 = _mm_loadu_ps (&(P(ix1, iy3, iz0)));  
  // cP
  tmp[0]   = _mm_mul_ps (cv, reg0);
  tmp[1]   = _mm_mul_ps (cv, reg1);
  tmp[2]   = _mm_mul_ps (cv, reg2);
  tmp[3]   = _mm_mul_ps (cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  cP[1]    = _mm_hadd_ps (tmp[0], tmp[1]);
  // dcP
  tmp[0]   = _mm_mul_ps (dcv, reg0);
  tmp[1]   = _mm_mul_ps (dcv, reg1);
  tmp[2]   = _mm_mul_ps (dcv, reg2);
  tmp[3]   = _mm_mul_ps (dcv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  dcP[1]   = _mm_hadd_ps (tmp[0], tmp[1]);
  // d2cP
  tmp[0]   = _mm_mul_ps (d2cv, reg0);
  tmp[1]   = _mm_mul_ps (d2cv, reg1);
  tmp[2]   = _mm_mul_ps (d2cv, reg2);
  tmp[3]   = _mm_mul_ps (d2cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  d2cP[1]  = _mm_hadd_ps (tmp[0], tmp[1]);

  reg0 = _mm_loadu_ps (&(P(ix2, iy0, iz0)));    
  reg1 = _mm_loadu_ps (&(P(ix2, iy1, iz0)));  
  reg2 = _mm_loadu_ps (&(P(ix2, iy2, iz0)));  
  reg3 = _mm_loadu_ps (&(P(ix2, iy3, iz0)));  
  // cP
  tmp[0]   = _mm_mul_ps (cv, reg0);
  tmp[1]   = _mm_mul_ps (cv, reg1);
  tmp[2]   = _mm_mul_ps (cv, reg2);
  tmp[3]   = _mm_mul_ps (cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  cP[2]    = _mm_hadd_ps (tmp[0], tmp[1]);
  // dcP
  tmp[0]   = _mm_mul_ps (dcv, reg0);
  tmp[1]   = _mm_mul_ps (dcv, reg1);
  tmp[2]   = _mm_mul_ps (dcv, reg2);
  tmp[3]   = _mm_mul_ps (dcv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  dcP[2]   = _mm_hadd_ps (tmp[0], tmp[1]);
  // d2cP
  tmp[0]   = _mm_mul_ps (d2cv, reg0);
  tmp[1]   = _mm_mul_ps (d2cv, reg1);
  tmp[2]   = _mm_mul_ps (d2cv, reg2);
  tmp[3]   = _mm_mul_ps (d2cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  d2cP[2]  = _mm_hadd_ps (tmp[0], tmp[1]);

  reg0 = _mm_loadu_ps (&(P(ix3, iy0, iz0)));    
  reg1 = _mm_loadu_ps (&(P(ix3, iy1, iz0)));  
  reg2 = _mm_loadu_ps (&(P(ix3, iy2, iz0)));  
  reg3 = _mm_loadu_ps (&(P(ix3, iy3, iz0)));  
  // cP
  tmp[0]   = _mm_mul_ps (cv, reg0);
  tmp[1]   = _mm_mul_ps (cv, reg1);
  tmp[2]   = _mm_mul_ps (cv, reg2);
  tmp[3]   = _mm_mul_ps (cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  cP[3]    = _mm_hadd_ps (tmp[0], tmp[1]);
  // dcP
  tmp[0]   = _mm_mul_ps (dcv, reg0);
  tmp[1]   = _mm_mul_ps (dcv, reg1);
  tmp[2]   = _mm_mul_ps (dcv, reg2);
  tmp[3]   = _mm_mul_ps (dcv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  dcP[3]   = _mm_hadd_ps (tmp[0], tmp[1]);
  // d2cP
  tmp[0]   = _mm_mul_ps (d2cv, reg0);
  tmp[1]   = _mm_mul_ps (d2cv, reg1);
  tmp[2]   = _mm_mul_ps (d2cv, reg2);
  tmp[3]   = _mm_mul_ps (d2cv, reg3);
  tmp[0]   = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1]   = _mm_hadd_ps (tmp[2], tmp[3]);
  d2cP[3]  = _mm_hadd_ps (tmp[0], tmp[1]);

  // Now compute bcP and bdcP
  tmp[0] = _mm_mul_ps (bv, cP[0]);
  tmp[1] = _mm_mul_ps (bv, cP[1]);
  tmp[2] = _mm_mul_ps (bv, cP[2]);
  tmp[3] = _mm_mul_ps (bv, cP[3]);

  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  bcP    = _mm_hadd_ps (tmp[0], tmp[1]);

  tmp[0] = _mm_mul_ps (bv, dcP[0]);
  tmp[1] = _mm_mul_ps (bv, dcP[1]);
  tmp[2] = _mm_mul_ps (bv, dcP[2]);
  tmp[3] = _mm_mul_ps (bv, dcP[3]);

  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  bdcP    = _mm_hadd_ps (tmp[0], tmp[1]);


  ///////////////////
  // Compute value //
  ///////////////////
  tmp[0] = _mm_mul_ps (av, bcP);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&val, tmp[0]);
  

  //////////////////////
  // Compute gradient //
  //////////////////////
  tmp[0] = _mm_mul_ps (dav, bcP);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(grad[0]), tmp[0]);
  
  tmp[0] = _mm_mul_ps (dbv, cP[0]);
  tmp[1] = _mm_mul_ps (dbv, cP[1]);
  tmp[2] = _mm_mul_ps (dbv, cP[2]);
  tmp[3] = _mm_mul_ps (dbv, cP[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[0] = _mm_mul_ps  (av, tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&grad[1], tmp[0]);
  
  tmp[0] = _mm_mul_ps (av, bdcP);
  tmp[0] = _mm_hadd_ps(tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps(tmp[0], tmp[0]);
  _mm_store_ss (&grad[2], tmp[0]);


  ////////////////////////////////
  // Compute second derivatives //
  ////////////////////////////////

  // (0,0) component
  tmp[0] = _mm_mul_ps (d2av, bcP);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(secDerivs(0,0)), tmp[0]);

  // (0,1) component
  tmp[0] = _mm_mul_ps (d2bv, cP[0]);
  tmp[1] = _mm_mul_ps (d2bv, cP[1]);
  tmp[2] = _mm_mul_ps (d2bv, cP[2]);
  tmp[3] = _mm_mul_ps (d2bv, cP[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[0] = _mm_mul_ps (dav, tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(secDerivs(0,1)), tmp[0]);
  secDerivs(1,0) = secDerivs(0,1);

  // (1,1) component
  tmp[0] = _mm_mul_ps (d2bv, cP[0]);
  tmp[1] = _mm_mul_ps (d2bv, cP[1]);
  tmp[2] = _mm_mul_ps (d2bv, cP[2]);
  tmp[3] = _mm_mul_ps (d2bv, cP[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[0] = _mm_mul_ps (av, tmp[0]);
  tmp[0] = _mm_mul_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_mul_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(secDerivs(1,1)), tmp[0]);

  // (1,2) component
  tmp[0] = _mm_mul_ps (dbv, dcP[0]);
  tmp[1] = _mm_mul_ps (dbv, dcP[1]);
  tmp[2] = _mm_mul_ps (dbv, dcP[2]);
  tmp[3] = _mm_mul_ps (dbv, dcP[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[0] = _mm_mul_ps (av, tmp[0]);
  tmp[0] = _mm_mul_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_mul_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(secDerivs(1,2)), tmp[0]);
  _mm_store_ss (&(secDerivs(2,1)), tmp[0]);

  // (2,2) component
  tmp[0] = _mm_mul_ps (bv, d2cP[0]);
  tmp[1] = _mm_mul_ps (bv, d2cP[1]);
  tmp[2] = _mm_mul_ps (bv, d2cP[2]);
  tmp[3] = _mm_mul_ps (bv, d2cP[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[1] = _mm_hadd_ps (tmp[2], tmp[3]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[1]);
  tmp[0] = _mm_mul_ps  (av, tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  tmp[0] = _mm_hadd_ps (tmp[0], tmp[0]);
  _mm_store_ss (&(secDerivs(2,2)), tmp[0]);
}

#endif




#endif
