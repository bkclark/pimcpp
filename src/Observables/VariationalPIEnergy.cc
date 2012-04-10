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

#include "VariationalPIEnergy.h"
#include <Common/MatrixOps/MatrixOps.h>



///Evalues the quantity <r|\exp(-\Beta*H)|r'>.  Will be fixed soon
///to evaluate correctly in periodic boundary conditions.
double 
VariationalPIEnergyClass::rho(int i, int j)
{
  double T=1.0/(PathData.Path.tau*PathData.Path.TotalNumSlices);
  double lambda=PathData.Path.Species(0).lambda;
  double total=0;
  for (int numX=-NumImages;numX<=NumImages;numX++){
    for (int numY=-NumImages;numY<=NumImages;numY++){
      for (int numZ=-NumImages;numZ<=NumImages;numZ++){
	dVec disp=PathData.Path(0,i)-PathData.Path(0,j);
	PathData.Path.PutInBox(disp);
	disp[0]=disp[0]+numX*PathData.Path.GetBox()[0];
	disp[1]=disp[1]+numY*PathData.Path.GetBox()[1];
	disp[2]=disp[2]+numZ*PathData.Path.GetBox()[2];
	double dist2=dot(disp,disp); //disp*disp;
	total+=exp(-T*dist2/(4.0*lambda));
      }
    }
  }
  return total;
}

  ///Returns the prefactor for the energy estimator
double VariationalPIEnergyClass::DRho(int i, int j)
{
  double total=0;
  double T=1.0/(PathData.Path.tau*PathData.Path.TotalNumSlices);
  double lambda=PathData.Path.Species(0).lambda;
  for (int numX=-NumImages;numX<=NumImages;numX++){
    for (int numY=-NumImages;numY<=NumImages;numY++){
      for (int numZ=-NumImages;numZ<=NumImages;numZ++){
	dVec disp=PathData.Path(0,i)-PathData.Path(0,j);
	PathData.Path.PutInBox(disp);
	disp[0]=disp[0]+numX*PathData.Path.GetBox()[0];
	disp[1]=disp[1]+numY*PathData.Path.GetBox()[1];
	disp[2]=disp[2]+numZ*PathData.Path.GetBox()[2];
	double dist2=dot(disp,disp); //disp*disp;
	//	total+=(3.0*T/2.0 - dist2*T*T/2)*exp(-T*dist2/2.0);
	total+=(3.0*T/2.0 - dist2*T*T/(4.0*lambda))*exp(-T*dist2/(4.0*lambda));
      }
    }
  }
  return total;
}



void VariationalPIEnergyClass::Accumulate()
{
  int NumParticles=PathData.Path.NumParticles();
  for (int ptcl1=0;ptcl1<NumParticles;ptcl1++){
    for (int ptcl2=0;ptcl2<NumParticles;ptcl2++){
      DetMatrix(ptcl1,ptcl2)=rho(ptcl1,ptcl2);
    }
  }
  double Rho = Determinant(DetMatrix);
  

  double TotalNum=0.0;
  double TotalDen=0.0;

  for (int n = 0; n<NumParticles;n++){
    for (int ptcl1=0;ptcl1<NumParticles;ptcl1++){
      for (int ptcl2=0;ptcl2<NumParticles;ptcl2++){
	if(ptcl1==n)
	  DetMatrix(ptcl1,ptcl2)=DRho(ptcl1,ptcl2);
	else
	  DetMatrix(ptcl1,ptcl2)=rho(ptcl1,ptcl2);
      }
    }
    double Det = Determinant(DetMatrix);
    TotalNum+= Det;
    TotalDen+= Rho;
  }
  Energy+= TotalNum/Rho;///TotalDen;
  //  EnergySq+= (TotalNum*TotalNum)/(Rho*Rho);
  
  NumSamples++;
  
}


void VariationalPIEnergyClass::WriteBlock()
{
  Energy=Energy/(double)NumSamples;
  EnergyVar.Write(Energy);
  EnergyVar.Flush();
  Energy=0.0;
  NumSamples = 0; 

}

void VariationalPIEnergyClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  int n=PathData.Path.NumParticles();
  DetMatrix.resize(n,n);
}
