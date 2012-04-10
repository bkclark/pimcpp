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

#include "../PathDataClass.h"
#include "ShortRangePotClass.h"

ShortRangePotClass::ShortRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

double 
ShortRangePotClass::V(int slice)
{
  double val = 0.0;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      double dist;
      dVec disp;
      
      PathData.Path.DistDisp(slice, ptcl1, ptcl2, dist, disp);
      PairActionFitClass& pa=
	*(PairMatrix(species1, PathData.Path.ParticleSpeciesNum(ptcl2)));
      
      val += pa.V(dist);
      if (pa.IsLongRange() && PathData.Actions.UseLongRange)
	val -= pa.Vlong(dist);
    }
  }
  return val;
}
