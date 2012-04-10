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

#include "SuperfluiDrop.h"

///////////////////////////////////////////////////////
//Superfluid estimator for a droplet = 2*m*<A(z)^2> ///
//                        /(beta*lambda*MomentInertia)///
//Vector A=0.5 Sum_{particle i,slice j}[rij x rij+1]///
//MomentInertia=<Sum{particle i slice j}[m*rij.rij+1]>///
///////////////////////////////////////////////////////




void SuperfluiDrop::WriteBlock()
{

  area=area/anorm;// Can I simply omit areasquared and replace it with area
  areaSquared.Write(area);//area=Write(area)
  areaSquared.Flush();
  area=0.0;
  anorm=0.0;

  momi=momi/mominorm;
  momInertia.Write(momi);
  momInertia.Flush();
  momi=0.0;
  mominorm=0.0;

  TotalCounts=0;

}


void SuperfluiDrop::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  string speciesName;
    assert(in.ReadVar("Species1",speciesName));
  for (int spec=0;spec<PathData.NumSpecies();spec++){ //???what is Species ?
    if (PathData.Species(spec).Name==speciesName){
      Species=spec;
    }
  }
  if (PathData.Path.Communicator.MyProc()==0){
    IOSection.WriteVar("Type","DropletSuperfluidity");
    IOSection.WriteVar("Cumulative", false);//???
  }
  area=0.0;
  momi=0.0;
  anorm=0.0;
  mominorm=0.0;
  TotalCounts=0;
  /// Now write the one-time output variables
//   if (PathData.Path.Communicator.MyProc()==0)
//     WriteInfo();

}

void SuperfluiDrop::Accumulate()
{
  TotalCounts++;
  PathClass &Path= PathData.Path;
  SpeciesClass &species=PathData.Path.Species(Species);

  dVec r;
  dVec angMom1,angMom2;
  dVec v;
  
  int nx,ny,nz,kjoinslice,permPtcl;
  double xold, yold, x, y;
  double a,moi;

  for (nz=0;nz<=2;nz++){// replace 3 with an input var
    nx=0;
    ny=1;

    if (0 == nz) nx=2;
    if (1 == nz) ny=2;

    a=0.0;
    moi=0.0;

    for (int ptcl=species.FirstPtcl;ptcl<=species.LastPtcl;ptcl++){

      kjoinslice=PathData.Join;//PathData.Path.Join;
      permPtcl=PathData.Path.Permutation(ptcl);

      xold=Path(kjoinslice,permPtcl)[nx];
      yold=Path(kjoinslice,permPtcl)[ny];

      for (int slice=0;slice<PathData.NumTimeSlices();slice++){
	kjoinslice++;
	if(kjoinslice == PathData.NumTimeSlices() ) 
	  kjoinslice=kjoinslice-PathData.NumTimeSlices();
	//	  cerr <<"kjoinslice "<<kjoinslice<<endl;


	x=Path( kjoinslice,ptcl)[nx];
	y=Path( kjoinslice,ptcl)[ny];

	a+= x*yold-y*xold;
	moi+=x*xold+y*yold ;

	xold=x;
	yold=y;
      }
    }
    a=0.5*a;
    area+= a*a;
    anorm+= 1.0;
    momi+= moi;
    mominorm+=PathData.NumTimeSlices();
  }
 cerr<<"PathData.NumTimeSlices(  "<<PathData.NumTimeSlices()<<endl;
  cerr<<"anorm "<<anorm<<" mominorm "<<mominorm<<endl;
  cerr <<" area " <<area<<" "<<area/anorm <<endl;
  cerr<<"MomI "<<moi<<" "<<momi/mominorm <<endl;
}
