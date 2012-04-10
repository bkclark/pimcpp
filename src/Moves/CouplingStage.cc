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

#include <Common/MPI/Communication.h>
#include "CouplingStage.h"
void CouplingStageClass::Read(IOSectionClass &in)
{

  LocalStageClass::Read(in);
}

void CouplingStageClass::WriteRatio()
{
  //  cerr<<"About to write my ratio"<<endl;
  Array<double,1> acceptRatio(2);
  acceptRatio(0)=(double)AcceptRatio(0);
  acceptRatio(1)=(double)AcceptRatio(1);
  AcceptRatioVar.Write(acceptRatio);
  AcceptRatioVar.Flush();
  //  cerr<<"done writing my ratio"<<endl;
}


void CouplingStageClass::Accept()
{
  //  cerr<<"I have accepted a lambda of "<<PathData.Path.ExistsCoupling<<endl;

}

void CouplingStageClass::Reject()
{
  //  cerr<<"I have rejected a lambda of "<<PathData.Path.ExistsCoupling<<endl;

}


///Chooses the time slices and moves the join so that the join is in
///the correct place for that time slice.
void CouplingStageClass::ChooseTimeSlices(int &slice1,int &slice2)
{
  //  slice1=0;
  //  slice2=PathData.NumTimeSlices();

}


// double CouplingStageClass::Sample(int &slice1,int &slice2, 
// 	      Array<int,1> &activeParticles)
// {
//   //slice1=0;
//   //    slice2=PathData.NumTimeSlices()-1;
//   //  activeParticles.resize(PathData.NumParticles());
//   //  for (int ptcl=0;ptcl<activeParticles.size();ptcl++)
//   //    activeParticles(ptcl)=ptcl;
//   //  if (PathData.Path.ExistsCoupling==1)
//   //    PathData.Path.ExistsCoupling=0;
//   //  else
//   //    PathData.Path.ExistsCoupling=1;
//   //  return 1.0;
  


//   ///For free energy
//   slice1=0;
//   slice2=PathData.NumTimeSlices()-1;
//     activeParticles.resize(1);
//     activeParticles(0)=PathData.Path.Species(1).FirstPtcl;
//     //    cerr<<"My lambda is currently"<<PathData.Path.ExistsCoupling<<endl;
//     //    int myProc=PathData.Path.Communicator.MyProc();
//     int myProc=PathData.InterComm.MyProc();
//     double lowVal=((double)(myProc))/100.0;
//     double highVal=((double)(myProc))/100.0+1.0/100.0;
//     //    cerr<<PathData.Path.Communicator.MyProc()<<" "<<"lowVal is "<<lowVal<<"and highVal is "<<highVal<<endl;
//     //    //    if (abs(PathData.Path.ExistsCoupling-lowVal)<1e-5){
//     //    //      PathData.Path.ExistsCoupling=highVal;
// //       //      cerr<<"Trying to switch from "<<PathData.Path.ExistsCoupling<<" to "<<highVal<<endl;
// //    }
// //    else{
// //      PathData.Path.ExistsCoupling=lowVal;
// //       //      cerr<<"Trying to switch from "<<PathData.Path.ExistsCoupling<<" to "<<lowVal<<endl;
// //    }
//     double oldCoupling=PathData.Path.ExistsCoupling;
//     if (PathData.Random.Local()>0.5)
//       PathData.Path.ExistsCoupling=PathData.Path.ExistsCoupling+PathData.Path.Random.Local()*0.1;
//     else
//       PathData.Path.ExistsCoupling=PathData.Path.ExistsCoupling-PathData.Path.Random.Local()*0.1;
//     if (PathData.Path.ExistsCoupling>1.0)
//       PathData.Path.ExistsCoupling=PathData.Path.ExistsCoupling-1.0;
//     if (PathData.Path.ExistsCoupling<0.0)
//       PathData.Path.ExistsCoupling=PathData.Path.ExistsCoupling+1.0;
//     //    cerr<<"My lambda is currently trying to be "<<PathData.Path.ExistsCoupling<<endl;
//     return oldCoupling/PathData.Path.ExistsCoupling;
//       //    return 1.0;
// }



double CouplingStageClass::Sample(int &slice1,int &slice2, 
	      Array<int,1> &activeParticles)
{
  int myProc=PathData.InterComm.MyProc();
  //  cerr<<"Inside coupling stage "<<myProc<<endl;
  ///For free energy
  slice1=0;
  slice2=PathData.NumTimeSlices()-1;
  activeParticles.resize(1);
  activeParticles(0)=PathData.Path.Species(1).FirstPtcl;
  double lowVal=((double)(myProc))/100.0;
  double highVal=((double)(myProc))/100.0+1.0/100.0;
  if (abs(PathData.Path.ExistsCoupling-lowVal)<1e-5)
    PathData.Path.ExistsCoupling=highVal;
  else
    PathData.Path.ExistsCoupling=lowVal;
  //  cerr<<"Leaving coupling stage "<<myProc<<endl;
  return 1.0;
}
