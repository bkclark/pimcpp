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

#include "PlaneDensity.h"

// ///Only works in cubic box
// int PlaneDensityClass::IntoGrid(double num)
// {
//   dVec box=PathData.Path.GetBox();
//   double maxLen=max(box[0],max(box[1],box[2]));
//   while (num>maxLen/2) 
//     num-=maxLen;
//   while (num<-maxLen/2)
//     num+=maxLen;
//   num+=maxLen/2.0;
//   int myNum=(int)floor(num/(maxLen/(double)(Grid.extent(0))));
//   //  //  cerr<<num;
//   //  cerr<<maxLen/(double)(Grid.extent(0));
//   //  cerr<<"My num is "<<myNum<<endl;
//   return myNum;
// }

int PlaneDensityClass::IntoGrid(double num,int dim)
{
  dVec box=PathData.Path.GetBox();
  double boxLen=box[dim];
  while (num>boxLen/2) 
    num-=boxLen;
  while (num<-boxLen/2)
    num+=boxLen;
  num+=boxLen/2.0;
  int myNum=(int)floor(num/(boxLen/(double)(Grid.extent(dim))));
  return myNum;
}


// Fix to include final link between link M and 0
void PlaneDensityClass::Accumulate()
{
  for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      NumSamples++;
      int nx=IntoGrid(PathData.Path(slice,ptcl)[0],0);
      int ny=IntoGrid(PathData.Path(slice,ptcl)[1],1);
      if (nx<NumSamples && ny<NumSamples)
	Grid(nx,ny)=Grid(nx,ny)+1;
      
    }
}

void PlaneDensityClass::WriteBlock()
{
  Array<double,2> SumGrid(Grid.extent(0),Grid.extent(1));
  //Path.Communicator.Sum(Grid, SumGrid);
  Path.Communicator.AllSum(Grid, SumGrid);
 
  double norm = 1.0/((double)NumSamples);
  SumGrid=SumGrid/norm;
  GridVar.Write(SumGrid);
  GridVar.Flush();
  NumSamples = 0;
  Grid=0.0;

}

void PlaneDensityClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);

  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Grid");
  }
  
}



void PlaneDensityClass::WriteInfo()
{


}
