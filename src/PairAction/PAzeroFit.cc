/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
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
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "PAzeroFit.h"

double PAzeroFitClass::U(double q, double z, double s2, int level)
{
  return 0.0;
}

double PAzeroFitClass::dU(double q, double z, double s2, int level)
{
  return 0.0;
}

double PAzeroFitClass::V(double r)
{
  return 0.0;
}

void 
PAzeroFitClass::Derivs (double q, double z, double s2, int level,
			double &d_dq, double &d_dz, double &d_ds)
{
  d_dq = 0.0;
  d_dz = 0.0;
  d_ds = 0.0;
}

bool PAzeroFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
  //  SmallestBeta = smallestBeta;
  // Read Particles;
    assert(in.OpenSection("Fits"));
    assert(in.OpenSection("Particle1"));
    Particle1.Read(in);
    in.CloseSection();
    assert(in.OpenSection("Particle2"));
    Particle2.Read(in);
    in.CloseSection();
  //  lambda = Particle1.lambda + Particle2.lambda;
  //  assert (lambda == 0.0);

  // Read Potential;
  //  assert(in.OpenSection("Potential"));
  //  Potential = ReadPH(in);
  //  in.CloseSection();
  //  in.CloseSection();
  return true;
}



bool PAzeroFitClass::IsLongRange()
{
  return false;
}
