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

#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"

class BisectionStageSphereClass : public LocalStageClass
{
public:
  void WriteRatio();
  ObservableDouble AcceptRatioVar;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  ///Really needs to go someplace else, but for the moment we will keep it here
  double SphereRadius; 
  void Accept();
  void Reject();
  void CartesianToSpherical(dVec &r, double &theta,double &phi);
  void SurfaceUnitVectors(dVec &Er,dVec &Etheta, dVec &Ephi, double theta, double phi);

  TinyVector<double,3> SphericalToCartesian(double &a,double &theta, double &phi);
  void ProjectOntoSphere(dVec &r, double a);
  void RotateAroundX(dVec &vec, double theta);
  void RotateAroundY(dVec &vec, double theta);
  void RotateAroundZ(dVec &vec, double theta);
  double angleInBetween(dVec &r1, dVec &r2);
  void Read(IOSectionClass &in);
  BisectionStageSphereClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection),
    AcceptRatioVar("AcceptRatio",OutSection,pathData.Path.Communicator) 
  { 
    //do nothing for now
    BisectionLevel = level;
    SphereRadius=2.0;
  }
};

#endif
