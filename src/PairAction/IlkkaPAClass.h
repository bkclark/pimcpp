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

#ifndef IlkkaPAClass_H
#define IlkkaPAClass_H

#include "PAFitBase.h"

///This is the pair action class. It uses the following formula in
///order to calculate the pair action
/*! \f[\frac{u(q,0;\tau)+\sum_{k=1}^n 
 * A_{j}(q;\tau) s^{2k}\f]   */
class IlkkaPAClass : public PairActionFitClass
{
 public:
  // Parameters
  string type1, type2;
  double Z1, Z2, tau;
  bool longRange;
  int nOrder;

  // Read values
  Array<double, 1> r_u, k_u, x_u, y_u, r_du, k_du, x_du, y_du, r_v, k_v;
  Array<double, 1> uShort_r, uLong_k, duShort_r, duLong_k, vShort_r, vLong_k;
  Array<double, 2> uOffDiag_xy, duOffDiag_xy;
  Array<Array<double, 1>, 1> r_A_u, A_u, r_A_du, A_du;
  Array<double,1> Potential;

  // Grids
  GeneralGrid r_u_grid, x_u_grid, y_u_grid;
  GeneralGrid r_du_grid, x_du_grid, y_du_grid;
  GeneralGrid r_v_grid;
  Array<GeneralGrid, 1> r_A_u_grid, r_A_du_grid;

  // Splines
  CubicSpline uShort_r_spline, duShort_r_spline, vShort_r_spline;
  Array<CubicSpline, 1> A_u_spline, A_du_spline;
  BicubicSpline uOffDiag_xy_spline, duOffDiag_xy_spline;

  // Constants
  double uLong_r0, uLong_k0, duLong_r0, duLong_k0, vLong_r0, vLong_k0;
  double kCutoff;

  /// Function to read Ilkka's squarer file input.
  inline bool Read(IOSectionClass &IOSection, double desiredTau, int numLevels);
  void ReadIlkkaHDF5(string fileName);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  double V(double r);
  bool IsLongRange();
};

inline bool IlkkaPAClass::Read(IOSectionClass &in, double x, int y)
{
  Name = "IlkkaPAClass";
  string fileName;
  assert(in.OpenSection("Fits"));
    assert(in.ReadVar("NumOffDiagonalTerms",nOrder));
    assert(in.OpenSection("Particle1"));
      Particle1.Read(in);
    in.CloseSection();
    assert(in.OpenSection("Particle2"));
      Particle2.Read(in);
    in.CloseSection();
    lambda = Particle1.lambda + Particle2.lambda;
    assert(in.ReadVar("Ilkkadmfile",fileName));
    if(!in.ReadVar("longRange",longRange))
      longRange = false;
    if(!in.ReadVar("vLongRange",vLongRange))
      vLongRange = false;
    ReadIlkkaHDF5(fileName.c_str());
  in.CloseSection();
  return true;
}


#endif
