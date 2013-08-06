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

#include "SpecificHeat.h"

void SpecificHeatClass::Accumulate()
{
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  NumSamples++;
  
  PairActionFitClass *specificHeat1;
  PairActionFitClass *specificHeat2;
  double tauValue1;
  double tauValue2;
  ///Get pair array for specific heat 1
  for (int i=0;i<PathData.Actions.SpecificHeatPairArray.size();i++){
    if (PathData.Actions.SpecificHeatPairArray(i)->Particle1.Name=="SpecificHeatHe1"){
      specificHeat1=PathData.Actions.SpecificHeatPairArray(i);
      tauValue1=PathData.Actions.TauValues(i);
      assert(PathData.Actions.SpecificHeatPairArray(i)->Particle2.Name=="SpecificHeatHe1");
    }
    if (PathData.Actions.SpecificHeatPairArray(i)->Particle1.Name=="SpecificHeatHe2"){
      specificHeat2=PathData.Actions.SpecificHeatPairArray(i);
      tauValue2=PathData.Actions.TauValues(i);
      assert(PathData.Actions.SpecificHeatPairArray(i)->Particle2.Name=="SpecificHeatHe2");
    }
  }
  Array<int,1> changedParticles;
  changedParticles.resize(PathData.Path.NumParticles());
  for (int i=0;i<changedParticles.size();i++)
    changedParticles(i)=i;

  double srA1=
    PathData.Actions.ShortRange.SingleActionForcedPairAction(0, PathData.Path.NumTimeSlices()-1,
							   *specificHeat1);
  double srA2 =
    PathData.Actions.ShortRange.SingleActionForcedPairAction(0, PathData.Path.NumTimeSlices()-1,
							   *specificHeat2);
  

  double kA1= 
    PathData.Actions.Kinetic.SingleActionForcedTau  (0, PathData.Path.NumTimeSlices()-1,
						  changedParticles, 
						  0, tauValue1);
  double kA2 =
    PathData.Actions.Kinetic.SingleActionForcedTau (0, PathData.Path.NumTimeSlices()-1,
						    changedParticles, 
						    0, tauValue2);
    
  double kE1=PathData.Actions.Kinetic.d_dBetaForcedTau (0,PathData.Path.NumTimeSlices()-1,
							 0,tauValue1);
  
  double kE2=PathData.Actions.Kinetic.d_dBetaForcedTau (0,PathData.Path.NumTimeSlices()-1,
							 0,tauValue2);

  double srE1=PathData.Actions.ShortRange.d_dBetaForcedPairAction(0, PathData.Path.NumTimeSlices()-1,
								  *specificHeat1);

  double srE2=PathData.Actions.ShortRange.d_dBetaForcedPairAction (0, 
								   PathData.Path.NumTimeSlices()-1,
								   *specificHeat2);

  double E1=srE1+kE1;
  double E2=exp(-((srA2+kA2)-(srA1+kA1)))*(srE2+kE2);
  EA+=E1;
  EB+=E2;
  WeightA+=1;
  WeightB+=exp(-((srA2+kA2)-(srA1+kA1)));
  cerr<<"KE1, KE2, srA1, srA2 :"<<kE1<<" "<<kE2<<" "<<srE1<<" "<<srE2<<" "<<WeightB<<endl;
  TotalSum+=(E1-E2);
  //  PathData.Actions.SpecificHeat (kinetic, dUShort, dULong, node, vShort, vLong);
  //PathData.Actions.SpecificHeat(Energies);
	//kinetic = Energies["kinetic"];
	//dUShort = Energies["dUShort"];
	//dULong = Energies["dULong"];
	//node = Energies["node"];
	//vShort = Energies["vShort"];
	//vLong = Energies["vLong"];

//   TotalSum   += kinetic + dUShort + dULong + node;// + tip5p;
//   KineticSum += kinetic;
//   dUShortSum += dUShort;
//   dULongSum  += dULong;
//   NodeSum    += node;
//   VShortSum  += vShort;
//   VLongSum   += vLong;

}

void SpecificHeatClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void SpecificHeatClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples);
  WeightAVar.Write(Prefactor*PathData.Path.Communicator.Sum(WeightA)*norm);
  WeightBVar.Write(Prefactor*PathData.Path.Communicator.Sum(WeightB)*norm);
  EAVar.Write(Prefactor*PathData.Path.Communicator.Sum(EA)*norm);
  EBVar.Write(Prefactor*PathData.Path.Communicator.Sum(EB)*norm);
  //  TotalVar.Write   (Prefactor*PathData.Path.Communicator.Sum(TotalSum)*norm);
  //  KineticVar.Write (Prefactor*PathData.Path.Communicator.Sum(KineticSum)*norm);
  //  dUShortVar.Write (Prefactor*PathData.Path.Communicator.Sum(dUShortSum)*norm);

  //  dULongVar.Write  (Prefactor*PathData.Path.Communicator.Sum(dULongSum)*norm);
  //  NodeVar.Write    (Prefactor*PathData.Path.Communicator.Sum(NodeSum)*norm);
  //  VShortVar.Write  (Prefactor*PathData.Path.Communicator.Sum(VShortSum)*norm);
  //  VLongVar.Write   (Prefactor*PathData.Path.Communicator.Sum(VLongSum)*norm);
  
  TotalSum       = 0.0;
  KineticSum     = 0.0;
  dUShortSum     = 0.0;
  dULongSum      = 0.0;
  NodeSum        = 0.0;
  VShortSum      = 0.0;
  VLongSum       = 0.0;
  EA=0.0;
  WeightA=0.0;
  EB=0.0;
  WeightB=0.0;
  NumSamples = 0;
}

void SpecificHeatClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}

