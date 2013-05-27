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

#include "QuinticSpline.h"
// #include "../config.h"

//HACK! // #define F77_QUINAT F77_FUNC(quinat,QUINAT)
//HACK! // extern "C" void 
//HACK! // F77_QUINAT (int *N, double X[], double Y[], 
//HACK! // 	    double B[], double C[], double D[], 
//HACK! // 	    double E[], double F[]);


// bool myIsNAN(double x)
// {
//   union doublechar {
//     double d;
//     unsigned char bytes[8];
//   } y;
//   y.d = x;
//   unsigned int exponent = ((y.bytes[1]&0x0f)<<8) | (y.bytes[0] & 0xfe);
//   if (exponent == 2047) {
//     unsigned long f0 = ((y.bytes[1]&0xf0)>>4) | (y.bytes[2] &0x0f);
//     unsigned long f1 = ((y.bytes[2]&0xf0)>>4) | (y.bytes[3] &0x0f);
//     unsigned long f2 = ((y.bytes[3]&0xf0)>>4) | (y.bytes[4] &0x0f);
//     unsigned long f3 = ((y.bytes[4]&0xf0)>>4) | (y.bytes[5] &0x0f);
//     unsigned long f4 = ((y.bytes[5]&0xf0)>>4) | (y.bytes[6] &0x0f);
//     unsigned long f5 = ((y.bytes[6]&0xf0)>>4) | (y.bytes[7] &0x0f);
//     unsigned long f6 = ((y.bytes[7]&0xf0)>>4);
//     unsigned long F = f0 | (f1<<8) | (f2<<16) | (f3<<24) | (f4<<32)
//       | (f5<<40) | (f6<<48);
//     if (F != 0)
//       return true;
//     else
//       return false;
//   }
//   else
//     return false;
// }


bool myIsNormal (double x)
{
  char out[100];
  snprintf (out, 50, "%f", x);
  bool notNormal = ((out[0]=='n') || (out[0]=='N')||
		    (out[0]=='i') || (out[1]=='i')||
		    (out[0]=='I') || (out[1]=='I'));
  return (!notNormal);
}



	  


void QuinticSpline::Update()
{
  //assert(1==2);
  ///QUINTIC SPLINE QUINAT NOT THERE!
  /// First, use double and triple knots to specify first and second
  /// derivatives at the boundary if we so desire.
  offset=0;
  FX(0) = (*grid)(0);
  FY(0) = Y(0);
//   cerr << "isnan(StartDeriv) = " << isnan(StartDeriv) << endl;
//   cerr << "isnnormal(StartDeriv) = " << isnormal(StartDeriv) << endl;
//   cerr << "StartDeriv        = " << StartDeriv << endl;
  if (!myIsNAN(StartDeriv)){
    offset++;
    FX(1) = (*grid)(0);
    FY(1) = StartDeriv;
    if (!myIsNAN(StartDeriv2)) {
      offset++;
      FX(2) = (*grid)(0);
      FY(2) = StartDeriv2;
    }
  }
  
  int i = grid->NumPoints + offset;
  if (!myIsNAN(EndDeriv)) {
    FX(i) = (*grid)(grid->NumPoints-1);
    FY(i) = EndDeriv;
    i++;
    if (!myIsNAN(EndDeriv2)) {
      FX(i) = (*grid)(grid->NumPoints-1);
      FY(i) = EndDeriv2;
    }
  }

  // Now fill in the rest of the values.
  for (int i=1; i<grid->NumPoints; i++) {
    FX(i+offset) = (*grid)(i);
    FY(i+offset) = Y(i);
  }
  int Fpoints = FX.size();
  // Call FORTRAN routine
  //HACK!  F77_QUINAT (&Fpoints, FX.data(), FY.data(), FB.data(), FC.data(),
  //HACK!		FD.data(), FE.data(), FF.data());

  // Now copy the data into our coefficents
  B(0) = FB(offset);
  C(0) = FC(offset);
  D(0) = FD(offset);
  E(0) = FE(offset);
  F(0) = FF(offset);
  for (int i=1; i<grid->NumPoints; i++) {
    B(i) = FB(i+offset);
    C(i) = FC(i+offset);
    D(i) = FD(i+offset);
    E(i) = FE(i+offset);
  }
  /// Last FF(N-1) is not set by quinat, so we don't copy to avoid a
  /// valgrind error.
  for (int i=1; i<(grid->NumPoints-1); i++)
    F(i) = FF(i+offset);
  F(grid->NumPoints-1) = 0.0;
  UpToDate=true;
}


void QuinticSpline::Write(IOSectionClass &outSection)
{
  outSection.WriteVar("StartDeriv", StartDeriv);
  outSection.WriteVar("EndDeriv", EndDeriv);
  outSection.WriteVar("StartDeriv2", StartDeriv2);
  outSection.WriteVar("EndDeriv2", EndDeriv2);
  outSection.WriteVar("Y", Y);
  
  outSection.NewSection("Grid");
  grid->Write(outSection);
  outSection.CloseSection();
}

void QuinticSpline::Read(IOSectionClass &inSection)
{
  assert(inSection.ReadVar("StartDeriv", StartDeriv));
  assert(inSection.ReadVar("EndDeriv", EndDeriv));
  assert(inSection.ReadVar("StartDeriv2", StartDeriv2));
  assert(inSection.ReadVar("EndDeriv2", EndDeriv2));
  Array<double,1> newY;
  assert(inSection.ReadVar("Y", newY));
  assert(inSection.OpenSection("Grid"));
  Grid *newGrid = ReadGrid(inSection);
  inSection.CloseSection();
  Init (newGrid, newY, StartDeriv, EndDeriv, StartDeriv2, EndDeriv2);
  Update();
}

QuinticSpline&
QuinticSpline::operator=(const QuinticSpline& spline)
{
  UpToDate = spline.UpToDate;
  Y.resize(spline.Y.shape());   Y  = spline.Y;
  FX.resize(spline.FX.shape()); FX = spline.FX;
  FY.resize(spline.FY.shape()); FY = spline.FY;
  FB.resize(spline.FB.shape()); FB = spline.FB;
  FC.resize(spline.FC.shape()); FC = spline.FC;
  FD.resize(spline.FD.shape()); FD = spline.FD;
  FE.resize(spline.FE.shape()); FE = spline.FE;
  FF.resize(spline.FF.shape()); FF = spline.FF;
  offset = spline.offset;
  StartDeriv = spline.StartDeriv;
  StartDeriv2 = spline.StartDeriv2;
  EndDeriv    = spline.EndDeriv;
  EndDeriv2   = spline.EndDeriv2;
  B.resize (spline.B.shape()); B = spline.B;
  C.resize (spline.C.shape()); C = spline.C;
  D.resize (spline.D.shape()); D = spline.D;
  E.resize (spline.E.shape()); E = spline.E;
  F.resize (spline.F.shape()); F = spline.F;
  I = spline.I;
  J = spline.J;
  grid = spline.grid;
  NumParams = spline.NumParams;

  return *this;
}
