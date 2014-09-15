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
  double rMax = uLong_r_spline.grid->End;
  if (r > rMax)
    r = rMax;
  if (rp > rMax)
    rp = rMax;
  double rMin = uLong_r_spline.grid->Start;
  if (r < rMin)
    r = rMin;
  if (rp < rMin)
    rp = rMin;

  double tmpU = u_xy_spline(x,y);

  // Add in off-diagonal part
  if (longRange)
    tmpU -= 0.5*(uLong_r_spline(r) + uLong_r_spline(rp));

  return tmpU;
}


double IlkkaPAClass::V(double r)
{
  // Limits
  double rMax = v_r_spline.grid->End;
  if (r >= rMax)
    return 0.;
  double rMin = v_r_spline.grid->Start;
  if (r < rMin)
    r = rMin;

  double tmpV = v_r_spline(r);

  // Add in off-diagonal part
  if (longRange)
    tmpV -= vLong_r_spline(r);

  return tmpV;
}


double IlkkaPAClass::dU(double q, double z, double s2, int level)
{
  double s = sqrt(s2);
  double x = q + 0.5*s;
  double y = q - 0.5*s;
  double r = q + 0.5*z;
  double rp = q - 0.5*z;

  // Limits
  double rMax = duLong_r_spline.grid->End;
  if (r > rMax)
    r = rMax;
  if (rp > rMax)
    rp = rMax;
  double rMin = duLong_r_spline.grid->Start;
  if (r < rMin)
    r = rMin;
  if (rp < rMin)
    rp = rMin;

  double tmpDU = du_xy_spline(x,y);

  // Add in off-diagonal part
  if (longRange)
    tmpDU -= 0.5*(duLong_r_spline(r) + duLong_r_spline(rp));

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
  if (longRange) {
    assert(h5In.OpenSection("diag"));
    assert(h5In.ReadVar("r",r_u));
    assert(h5In.ReadVar("u_r",u_r));
    assert(h5In.ReadVar("rLong",r_uLong));
    assert(h5In.ReadVar("uLong_r",uLong_r));
    assert(h5In.ReadVar("uLong_r0",uLong_r0));
    assert(h5In.ReadVar("k",k_u));
    assert(h5In.ReadVar("uLong_k",uLong_k));
    assert(h5In.ReadVar("uLong_k0",uLong_k0));
    h5In.CloseSection();
  }
  assert(h5In.OpenSection("offDiag"));
  assert(h5In.ReadVar("x",x_u));
  assert(h5In.ReadVar("y",y_u));
  assert(h5In.ReadVar("uOffDiag",uOffDiag_xy));
  assert(h5In.ReadVar("u_xy",u_xy));
  h5In.CloseSection();
  h5In.CloseSection();

  // Spline u
  if(longRange) {
    r_u_grid.Init(r_u);
    u_r_spline.Init(&r_u_grid, u_r);
    r_uLong_grid.Init(r_uLong);
    uLong_r_spline.Init(&r_uLong_grid, uLong_r);
  }
  x_u_grid.Init(x_u);
  y_u_grid.Init(y_u);
  uOffDiag_xy_spline.Init(&x_u_grid, &y_u_grid, uOffDiag_xy);
  u_xy_spline.Init(&x_u_grid, &y_u_grid, u_xy);

  // Read in du
  if (longRange) {
    assert(h5In.OpenSection("du"));
    assert(h5In.OpenSection("diag"));
    assert(h5In.ReadVar("r",r_du));
    assert(h5In.ReadVar("du_r",du_r));
    assert(h5In.ReadVar("rLong",r_duLong));
    assert(h5In.ReadVar("duLong_r",duLong_r));
    assert(h5In.ReadVar("duLong_r0",duLong_r0));
    assert(h5In.ReadVar("k",k_du));
    assert(h5In.ReadVar("duLong_k",duLong_k));
    assert(h5In.ReadVar("duLong_k0",duLong_k0));
    h5In.CloseSection();
  }
  assert(h5In.OpenSection("offDiag"));
  assert(h5In.ReadVar("x",x_du));
  assert(h5In.ReadVar("y",y_du));
  assert(h5In.ReadVar("duOffDiag",duOffDiag_xy));
  assert(h5In.ReadVar("du_xy",du_xy));
  h5In.CloseSection();
  h5In.CloseSection();

  // Spline du
  if (longRange) {
    r_du_grid.Init(r_du);
    du_r_spline.Init(&r_du_grid, du_r);
    r_duLong_grid.Init(r_duLong);
    duLong_r_spline.Init(&r_duLong_grid, duLong_r);
  }
  x_du_grid.Init(x_du);
  y_du_grid.Init(y_du);
  duOffDiag_xy_spline.Init(&x_du_grid, &y_du_grid, duOffDiag_xy);
  du_xy_spline.Init(&x_du_grid, &y_du_grid, du_xy);

  // Read in v
  assert(h5In.OpenSection("v"));
  assert(h5In.OpenSection("diag"));
  assert(h5In.ReadVar("r",r_v));
  assert(h5In.ReadVar("v_r",v_r));
  if (longRange) {
    assert(h5In.ReadVar("rLong",r_vLong));
    assert(h5In.ReadVar("vLong_r",vLong_r));
    assert(h5In.ReadVar("vLong_r0",vLong_r0));
    assert(h5In.ReadVar("k",k_v));
    assert(h5In.ReadVar("vLong_k",vLong_k));
    assert(h5In.ReadVar("vLong_k0",vLong_k0));
  }
  h5In.CloseSection();
  h5In.CloseSection();

  // Spline v
  r_v_grid.Init(r_v);
  v_r_spline.Init(&r_v_grid, v_r);
  if (longRange) {
    r_vLong_grid.Init(r_vLong);
    vLong_r_spline.Init(&r_vLong_grid, vLong_r);
  }

}
