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

#include "IlkkaPAClass.h"


double IlkkaPAClass::U (double q, double z, double s2, int level)
{
  double s = sqrt(s2);
  double x = q + 0.5*s;
  double y = q - 0.5*s;
  double r = q + 0.5*z;
  double rp = q - 0.5*z;

  // Limits
  double rMax = uShort_r_spline.grid->End;
  if (x > rMax)
    x = rMax;
  if (y > rMax)
    y = rMax;
  if (r > rMax)
    r = rMax;
  if (rp > rMax)
    rp = rMax;
  double rMin = uShort_r_spline.grid->Start;
  if (q < rMin)
    q = rMin;
  if (x < rMin)
    x = rMin;
  if (y < rMin)
    y = rMin;
  if (r < rMin)
    r = rMin;
  if (rp < rMin)
    rp = rMin;

  double tmpU = 0.;
  if (q <= rMax) {
    // Start with end-point action
    tmpU += 0.5*(uShort_r_spline(r) + uShort_r_spline(rp));
    //tmpU += uShort_r_spline(q);
  }

  // Add in off-diagonal part
  if (nOrder == -1) {
    tmpU += u_xy_spline(x,y);
    if (longRange)
      tmpU -= 0.5*(u_r_spline(r) + u_r_spline(rp));
  } else
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder)
      tmpU += A_u_spline(iOrder-1)(q) * pow(s2,iOrder);

  // Subtract off potential
  if (vLongRange) {
    rMin = vShort_r_spline.grid->Start;
    if (r < rMin)
      r = rMin;
    if (rp < rMin)
      rp = rMin;
    double tmpV = 0.5*Z1*Z2*tau*((1./r) + (1./rp));
    //double tmpV = 0.5*tau*(vShort_r_spline(r) + vShort_r_spline(rp));
    tmpU -= tmpV;
  }

  return tmpU;
}


double IlkkaPAClass::V(double r)
{
  // Limits
  double rMax = vShort_r_spline.grid->End;
  if (r >= rMax)
    return 0.;
  double rMin = vShort_r_spline.grid->Start;
  if (r < rMin)
    r = rMin;

  return vShort_r_spline(r);
}


double IlkkaPAClass::dU(double q, double z, double s2, int level)
{
  double s = sqrt(s2);
  double x = q + 0.5*s;
  double y = q - 0.5*s;
  double r = q + 0.5*z;
  double rp = q - 0.5*z;

  // Limits
  double rMax = duShort_r_spline.grid->End;
  if (x > rMax)
    x = rMax;
  if (y > rMax)
    y = rMax;
  if (r > rMax)
    r = rMax;
  if (rp > rMax)
    rp = rMax;
  double rMin = duShort_r_spline.grid->Start;
  if (q < rMin)
    q = rMin;
  if (x < rMin)
    x = rMin;
  if (y < rMin)
    y = rMin;
  if (r < rMin)
    r = rMin;
  if (rp < rMin)
    rp = rMin;

  double tmpDU = 0.;
  if (q <= rMax) {
    // Start with end-point action
    tmpDU += 0.5*(duShort_r_spline(r) + duShort_r_spline(rp));
    //tmpDU += duShort_r_spline(q);
  }

  // Add in off-diagonal part
  if (nOrder == -1) {
    tmpDU += du_xy_spline(x,y);
    if (longRange)
      tmpDU -= 0.5*(du_r_spline(r) + du_r_spline(rp));
  } else
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder)
      tmpDU += A_du_spline(iOrder-1)(q) * pow(s2,iOrder);

  // Subtract off potential
  if (vLongRange) {
    rMin = vShort_r_spline.grid->Start;
    if (r < rMin)
      r = rMin;
    if (rp < rMin)
      rp = rMin;
    double tmpV = 0.5*Z1*Z2*((1./r) + (1./rp));
    tmpDU -= tmpV;
  }

  return tmpDU;
}


bool IlkkaPAClass::IsLongRange()
{
  return false;
}


void IlkkaPAClass::ReadIlkkaHDF5(string fileName)
{
  // Read in hdf5 file
  IOSectionClass h5In;
  if(!h5In.OpenFile(fileName)) {
    cerr << "ERROR: Could not find pair action file " << fileName << ". Aborting..." << endl;
    abort();
  }

  // Read in info
  assert(h5In.OpenSection("Info"));
  assert(h5In.ReadVar("Z1",Z1));
  assert(h5In.ReadVar("Z2",Z2));
  assert(h5In.ReadVar("tau",tau));
  h5In.CloseSection();

  // Read in u
  assert(h5In.OpenSection("u"));
  assert(h5In.OpenSection("diag"));
  assert(h5In.ReadVar("r",r_u));
  assert(h5In.ReadVar("uShort_r",uShort_r));
  assert(h5In.ReadVar("u_r",u_r));
  if (longRange) {
    assert(h5In.ReadVar("uLong_r0",uLong_r0));
    assert(h5In.ReadVar("k",k_u));
    assert(h5In.ReadVar("uLong_k",uLong_k));
    assert(h5In.ReadVar("uLong_k0",uLong_k0));
  }
  h5In.CloseSection();
  assert(h5In.OpenSection("offDiag"));
  assert(h5In.ReadVar("x",x_u));
  assert(h5In.ReadVar("y",y_u));
  assert(h5In.ReadVar("uOffDiag",uOffDiag_xy));
  assert(h5In.ReadVar("u_xy",u_xy));
  if (nOrder > 0) {
    r_A_u.resize(nOrder);
    A_u.resize(nOrder);
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder) {
      ostringstream stream;
      stream << "A." << iOrder;
      string SectionTitle(stream.str());
      assert(h5In.OpenSection(SectionTitle));
      assert(h5In.ReadVar("r",r_A_u(iOrder-1)));
      assert(h5In.ReadVar("A",A_u(iOrder-1)));
      h5In.CloseSection();
    }
  }
  h5In.CloseSection();
  h5In.CloseSection();

  // Spline u
  r_u_grid.Init(r_u);
  uShort_r_spline.Init(&r_u_grid, uShort_r);
  x_u_grid.Init(x_u);
  y_u_grid.Init(y_u);
  u_r_spline.Init(&x_u_grid, u_r);
  uOffDiag_xy_spline.Init(&x_u_grid, &y_u_grid, uOffDiag_xy);
  u_xy_spline.Init(&x_u_grid, &y_u_grid, u_xy);
  if (nOrder > 0) {
    r_A_u_grid.resize(nOrder);
    A_u_spline.resize(nOrder);
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder) {
      r_A_u_grid(iOrder-1).Init(r_A_u(iOrder-1));
      A_u_spline(iOrder-1).Init(&r_A_u_grid(iOrder-1), A_u(iOrder-1));
    }
  }

  // Read in du
  assert(h5In.OpenSection("du"));
  assert(h5In.OpenSection("diag"));
  assert(h5In.ReadVar("r",r_du));
  assert(h5In.ReadVar("duShort_r",duShort_r));
  assert(h5In.ReadVar("du_r",du_r));
  if (longRange) {
    assert(h5In.ReadVar("duLong_r0",duLong_r0));
    assert(h5In.ReadVar("k",k_du));
    assert(h5In.ReadVar("duLong_k",duLong_k));
    assert(h5In.ReadVar("duLong_k0",duLong_k0));
  }
  h5In.CloseSection();
  assert(h5In.OpenSection("offDiag"));
  assert(h5In.ReadVar("x",x_du));
  assert(h5In.ReadVar("y",y_du));
  assert(h5In.ReadVar("duOffDiag",duOffDiag_xy));
  assert(h5In.ReadVar("du_xy",du_xy));
  if (nOrder > 0) {
    r_A_du.resize(nOrder);
    A_du.resize(nOrder);
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder) {
      ostringstream stream;
      stream << "A." << iOrder;
      string SectionTitle(stream.str());
      assert(h5In.OpenSection(SectionTitle));
      assert(h5In.ReadVar("r",r_A_du(iOrder-1)));
      assert(h5In.ReadVar("A",A_du(iOrder-1)));
      h5In.CloseSection();
    }
  }
  h5In.CloseSection();
  h5In.CloseSection();

  // Spline du
  r_du_grid.Init(r_du);
  duShort_r_spline.Init(&r_du_grid, duShort_r);
  x_du_grid.Init(x_du);
  y_du_grid.Init(y_du);
  du_r_spline.Init(&x_du_grid, du_r);
  duOffDiag_xy_spline.Init(&x_du_grid, &y_du_grid, duOffDiag_xy);
  du_xy_spline.Init(&x_du_grid, &y_du_grid, du_xy);
  if (nOrder > 0) {
    r_A_du_grid.resize(nOrder);
    A_du_spline.resize(nOrder);
    for (int iOrder=1; iOrder<nOrder+1; ++iOrder) {
      r_A_du_grid(iOrder-1).Init(r_A_du(iOrder-1));
      A_du_spline(iOrder-1).Init(&r_A_du_grid(iOrder-1), A_du(iOrder-1));
    }
  }

  // Read in v
  assert(h5In.OpenSection("v"));
  assert(h5In.OpenSection("diag"));
  assert(h5In.ReadVar("r",r_v));
  assert(h5In.ReadVar("vShort_r",vShort_r));
  if (longRange||vLongRange) {
    assert(h5In.ReadVar("vLong_r0",vLong_r0));
    assert(h5In.ReadVar("k",k_v));
    assert(h5In.ReadVar("vLong_k",vLong_k));
    assert(h5In.ReadVar("vLong_k0",vLong_k0));
  }
  h5In.CloseSection();
  h5In.CloseSection();

  // If using alternative long range, set u and du
  if (vLongRange) {
    uLong_r0 = tau*vLong_r0;
    duLong_r0 = vLong_r0;
    k_u.resize(k_v.size());
    k_du.resize(k_v.size());
    uLong_k.resize(k_v.size());
    duLong_k.resize(k_v.size());
    for (int i=0; i<k_v.size(); i++) {
      k_u(i) = k_v(i);
      k_du(i) = k_v(i);
      uLong_k(i) = tau*vLong_k(i);
      duLong_k(i) = vLong_k(i);
    }
    uLong_k0 = tau*vLong_k0;
    duLong_k0 = vLong_k0;
  }

  // Spline v
  r_v_grid.Init(r_v);
  vShort_r_spline.Init(&r_v_grid, vShort_r);

}
