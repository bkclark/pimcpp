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

#include "../PathDataClass.h"
#include "Tether.h"

std::string 
TetherClass::GetName()
{
  return "Tether";
}

///This has to be called after pathdata knows how many
///particles it has
void TetherClass::Read(IOSectionClass& in)
{
  string speciesString;
  assert(in.ReadVar("Species",speciesString));
  SpeciesNum=PathData.Path.SpeciesNum(speciesString);
  assert(SpeciesNum!=-1);
  cerr<<"My species num is "<<SpeciesNum<<endl;
  cerr<<"Speices num of he is "<<PathData.Path.SpeciesNum("He");
  assert(in.ReadVar("TetherCutoff",TetherCutoff));
  Array<double,2> positions;
  assert(in.ReadVar("TetheredSites",positions));
  TetheredSites.resize(PathData.Path.Species(SpeciesNum).LastPtcl-PathData.Path.Species(SpeciesNum).FirstPtcl+1);
  ///Verify you used the right number of points to compare against
  assert(positions.extent(0)==TetheredSites.size());
  assert(positions.extent(1)==NDIM);
  dVec pos;
  for (int loc=0;loc<TetheredSites.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    TetheredSites(loc) = pos;
  }      
  

}

TetherClass::TetherClass(PathDataClass &pathData ) : 
   ActionBaseClass (pathData)
 {
 }

double 
TetherClass::SingleAction (int slice1, int slice2,
				       const Array<int,1> &changedParticles, 
				       int level)
{
  int firstPtcl=PathData.Path.Species(SpeciesNum).FirstPtcl;
  //  cerr<<"My first ptcl is "<<firstPtcl<<endl;
  for (int slice=slice1;slice<=slice2;slice++){
    for (int ptclIndex=0;ptclIndex<changedParticles.size();ptclIndex++){
      int ptcl=changedParticles(ptclIndex);
      dVec disp=TetheredSites(ptcl-firstPtcl)-PathData.Path(slice,ptcl);
      PathData.Path.PutInBox(disp);
      double dist=sqrt(dot(disp,disp));
      if (dist>TetherCutoff)
	return 999999999.0;
    }
  }
  return 0.0;
}


double TetherClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}
