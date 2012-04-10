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

//#include "Common/Blitz.h"
#include "../PathDataClass.h"
#ifndef RANDOM_PERM_CLASS_H
#define RANDOM_PERM_CLASS_H


class RandomPermClass
{
 private:
/*   double GetActionEstimate(int i, int j, double fourLambdaBetaInv, */
/* 			   Array<double,2> speciesData); */
  Array<double,2> ActionPairs;




 public:
  void GetActionPairs();
  void Apply();
  void Read(IOSectionClass &inSection);
  PathDataClass &PathData;
  RandomPermClass(PathDataClass &myPathData);
  int Slice1,Slice2;
  int SpeciesNum;
  Array<int,1> PermArray; //Must resize perm array
  void GeneratePerm();
  Array<int,1> CurrentParticles;
};




/* inline double RandomPermClass::GetActionEstimate(int i,int j, */
/* 						 double fourLambdaBetaInv, */
/* 						 Array<double,2> speciesData) */
/* { */
/*   dVec disp_ij=speciesData(Slice2,j)-speciesData(Slice1,i); */
/*   PathData.Path.PutInBox(disp_ij); */
/*   double dist_ij=dot(disp_ij,disp_ij); */
/*   return exp(-dist_ij*fourLambdaBetaInv); */
/* } */





#endif
