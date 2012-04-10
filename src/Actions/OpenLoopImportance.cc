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
#include "OpenLoopImportance.h"

///BUG: Will not currently work in parallel


///This has to be called after pathdata knows how many
///particles it has
void OpenLoopImportanceClass::Read(IOSectionClass& in)
{
  string impChoiceString;
  string XisString;
  if (PathData.Path.OpenPaths){
    assert(in.OpenSection("OpenLoop"));

    assert(in.ReadVar("Xis",XisString));
    if (XisString=="Displacement")
      Xis=DISP;
    else if (XisString=="MinImageDisplacement")
      Xis=MIN_IMAGE_DISP;
    else if (XisString=="Distance")
      Xis=DIST;
    else if (XisString=="MinImageDistance")
      Xis=MIN_IMAGE_DIST;
    assert(in.ReadVar("ImportanceSample",impChoiceString));
    cerr<<"The ImportanceSampling string was "<<impChoiceString<<endl;
    if (impChoiceString=="None")
      ImpChoice=NOIMP;
    else if (impChoiceString=="Distance")
      ImpChoice=DISTIMP;
    else if (impChoiceString=="Displacement")
      ImpChoice=DISPXIMP;
    else if (impChoiceString=="Tuned")
      ImpChoice=TUNEDFUNCTION;
    else if (impChoiceString=="ReTuned")
      ImpChoice=RETUNEDFUNCTION;
    else if (impChoiceString=="Polynomial"){
      assert(in.ReadVar("Polynomial",Polynom));
      ImpChoice=POLYNOMIAL;
    }
    else if (impChoiceString=="Exponential"){
      assert(in.ReadVar("a",a));
      assert(in.ReadVar("s",s));
      assert(in.ReadVar("alpha",alpha));
      ImpChoice=EXPONENTIAL;
    }
    else {
      cerr<<"You have given an invalid choice"<<endl;
      assert(1==2);
    }
//     if (ImpChoice!=NOIMP){
//       string shiftTypeString;
//       assert(in.ReadVar("ShiftType",shiftTypeString));
//       if (shiftTypeString=="FixedShift")
// 	assert(in.ReadVar("ShiftAmount",Shift));
//       else if (shiftTypeString=="ProcShift"){
// 	int myProc=PathData.GetCloneNum();
// 	Shift=(myProc%16)+0.5;
//       }
//       else {
// 	cerr<<"I don't know what shift you want me to use"<<endl;
// 	assert(1==2);
//       }
//     }
    in.CloseSection();
  }
  else{
    cerr<<"You have an OpenLoop Action but haven't turned on open loops. This is not legal. Please fix your input file.";
    assert(1==2);
  }
}

OpenLoopImportanceClass::OpenLoopImportanceClass(PathDataClass &pathData ) : 
   ActionBaseClass (pathData)
 {
 }

double 
OpenLoopImportanceClass::SingleAction (int slice1, int slice2,
				       const Array<int,1> &changedParticles, 
				       int level)
{
  //  cerr<<"Open loop action begins"<<endl;
   int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
   ///If you don't own the open slice then just return. This is to
   ///make things happier with parallelization but it's not clear that
   ///this really works with parallelization right this moment.  
   if (procWithRefSlice != PathData.Path.Communicator.MyProc()) {
     return 0.0;
   }

   int openLink=PathData.Path.OpenLink;
   int openPtcl=PathData.Path.OpenPtcl;
   double x;
   if (Xis==MIN_IMAGE_DISP){
     dVec disp;
     double dist;
     PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			    dist,disp); //This is distance between
					//head and tail!
     x=disp(0);
   }
   else if (Xis==MIN_IMAGE_DIST){
     dVec disp;
     double dist;
     PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			    dist,disp); //This is distance between
					//head and tail!
     x=dist;
   }
   else if (Xis==DISP){
     int numLinks=PathData.Path.NumTimeSlices()-1;
     dVec disp=0.0;
     int currSlice=openLink;
     int currPtcl=openPtcl;
     int nextSlice=-1;
     int nextPtcl=-1;
     dVec linkDisp;
     while (nextSlice!=openLink || nextPtcl!=openPtcl){
       nextSlice = (currSlice + 1) % PathData.Path.NumTimeSlices();
       if (currSlice==PathData.Join)
	 nextPtcl=PathData.Path.Permutation(currPtcl);
       else 
	 nextPtcl=currPtcl;
       linkDisp=PathData.Path.VelocityBetweenPtcl(currSlice,currPtcl,
						  nextSlice,nextPtcl);
       disp=disp+linkDisp;
       currSlice=nextSlice;
       currPtcl=nextPtcl;
     }
     double dist=sqrt(dot(disp,disp));
     x=disp(0);
   }
   else if (Xis==DIST){
     int numLinks=PathData.Path.NumTimeSlices()-1;
     dVec disp=0.0;
     int currSlice=openLink;
     int currPtcl=openPtcl;
     int nextSlice=-1;
     int nextPtcl=-1;
     dVec linkDisp;
     while (nextSlice!=openLink || nextPtcl!=openPtcl){
       nextSlice = (currSlice + 1) % PathData.Path.NumTimeSlices();
       if (currSlice==PathData.Join)
	 nextPtcl=PathData.Path.Permutation(currPtcl);
       else 
	 nextPtcl=currPtcl;
       linkDisp=PathData.Path.VelocityBetweenPtcl(currSlice,currPtcl,
						  nextSlice,nextPtcl);
       disp=disp+linkDisp;
       currSlice=nextSlice;
       currPtcl=nextPtcl;
     }
     double dist=sqrt(dot(disp,disp));
     x=dist;
   }
   else {
     cerr<<"You need to specify how to establish what x is inside of the open action loop. Now aborting"<<endl;
     assert(1==2);
   }
   if (ImpChoice==NOIMP){
     return 0.0;
   }
   if (ImpChoice==POLYNOMIAL){
     double polyVal=0.0;
     double xToPwr=1.0;
     for (int pwr=0;pwr<Polynom.size();pwr++){
       polyVal += xToPwr*Polynom(pwr);
       xToPwr*=x;
     }
     return polyVal;
   }
   ///Returns a*exp(-alpha*(x-s)^2)
   else if (ImpChoice==EXPONENTIAL){
     //     cerr<<"open loop action ending"<<endl;
     //     cerr<<"Action info: "<<a<<" "<<x<<" "<<s<<" "<<alpha<<endl;
     return -log(a*exp(-alpha*(x-s)*(x-s)));
   }
   else {
     cerr<<"You have chosen a non-implemented importance sampling. Please choose again"<<endl;
     assert(1==2);
   }
}


// double 
// OpenLoopImportanceClass::SingleAction (int slice1, int slice2,
// 				       const Array<int,1> &changedParticles, 
// 				       int level)
// {

//   //    if (PathData.Path.Communicator.MyProc()==1)
//   //    cerr<<"Into importance action"<<endl;


//    int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
//    if (procWithRefSlice != PathData.Path.Communicator.MyProc()) {
//      //    if (PathData.Path.Communicator.MyProc()==1)
//      //      cerr<<"leaving imp action with 0"<<endl;

//      return 0.0;
//    }

//    int openLink=PathData.Path.OpenLink;
//    int openPtcl=PathData.Path.OpenPtcl;
//    //  if (PathData.Path.Communicator.MyProc()==1){
//    //    cerr<<"not returning with 0"<<endl;
//    //    cerr<<"Proc with ref slice is "<<procWithRefSlice<<endl;
//    //    cerr<<"The open link is "<<openLink<<endl;
//    //    cerr<<"The open ptcl is "<<openPtcl<<endl;
//    //    cerr<<"The values are "<<PathData.Path(openLink,PathData.Path.NumParticles())<<" and "<<PathData.Path(openLink,openPtcl)<<endl;
//    //  }
//    dVec disp;
//    double dist;
//    PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
// 			  dist,disp); //This is distance between head and tail!
//    //  int myProc=(PathData.InterComm.MyProc() % 16);
//    //  int myProc=PathData.GetCloneNum();
//    //  double shift=(myProc%16)+0.5;
//    //  //  shift=6.5;
//    //    if (PathData.Path.Communicator.MyProc()==1)
//    //      cerr<<"about to start hte ifs"<<endl;

//    if (ImpChoice==NOIMP){ 
//      //    if (PathData.Path.Communicator.MyProc()==1)
//      //      cerr<<"Into none action"<<endl;

//      //cerr<<"I have chosen no importance function"<<endl;
//      //    assert(1==3);
//      return 0.0;
//      //    return -log(1.0/(dist*dist+0.05)); 
//    }
//    else if (ImpChoice==DISTIMP){
//      //    if (PathData.Path.Communicator.MyProc()==1)
//      //      cerr<<"Into dist action"<<endl;
//      //    cerr<<"I have chosen a distance importance function"<<endl;
//      return -log(exp(-(dist-Shift)*(dist-Shift)));
//    }
//    else if (ImpChoice==DISPXIMP){
//      //    cerr<<"I have chosen a displacement importance function"<<endl;
//      //    if (PathData.Path.Communicator.MyProc()==1)
//      //      cerr<<"Into dispximp action"<<endl;



//      //    int numLinks=PathData.Path.NumTimeSlices()-1;
//      //    disp=0.0;
//      //    for (int slice=0;slice<numLinks;slice++) {
//      //      int realSlice=(openLink+slice) % numLinks;
//      //      int realSlicep1=(openLink+slice+1) % numLinks;
//      //      dVec linkDisp;
//      //      linkDisp=PathData.Path.Velocity(realSlice,realSlicep1,openPtcl);
//      //      disp =disp+linkDisp;
//      //    }
//      //    double  dist=sqrt(dot(disp,disp));

//      dVec disp=0.0;
//      double dist2;
//      int openLink=(int)(PathData.Path.OpenLink);
//      int openPtcl=(int)(PathData.Path.OpenPtcl);
//      PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
// 			    dist2,disp); //This is distance between head and tail!



//    int numLinks=PathData.Path.NumTimeSlices()-1;
//    disp=0.0;
//    int currSlice=openLink;
//    int currPtcl=openPtcl;
//    int nextSlice=-1;
//    int nextPtcl=-1;


//    dVec linkDisp;
//    //  cerr<<"hello"<<endl;
//    while (nextSlice!=openLink || nextPtcl!=openPtcl){
//      //    cerr<<nextSlice<<" "<<currSlice<<" "<<openLink<<endl;
//      //    cerr<<nextPtcl<<" "<<currPtcl<<" "<<openPtcl<<endl;

//      nextSlice = (currSlice + 1) % PathData.Path.NumTimeSlices();
//      //    if (nextSlice==0)
//      //      nextSlice=numLinks+1;
//      if (currSlice==PathData.Join)
//        nextPtcl=PathData.Path.Permutation(currPtcl);
//      else 
//        nextPtcl=currPtcl;
//      linkDisp=PathData.Path.VelocityBetweenPtcl(currSlice,currPtcl,nextSlice,nextPtcl);
//      disp=disp+linkDisp;
//      currSlice=nextSlice;
//      currPtcl=nextPtcl;
//    }
//    //HACK!
//    //  disp=disp-linkDisp;
//    double dist=sqrt(dot(disp,disp));

//   //  if (((dist-dist2)/PathData.Path.GetBox()[0])-floor((dist-dist2)/PathData.Path.GetBox()[0]+0.1)<1e-5)
//   //  cerr<<"Open Loop: dist, dist2, diff: "<<dist<<" "<<dist2<<" "<<dist-dist2<<endl;
//   //     return log(10.0)*(-0.29*dist*dist+0.0012*dist*dist*dist*dist); //-log(1.0/(dist*dist+0.05));
//   ////  if (dist<12)
//   ///    return 5*log(10.0)-2.75*dist;
//     //    return 5*log(10.0)-1.75*dist;
//   ////  else return 9999999;
//    //  double Shift=9.0;
//    //  double Shift=5.5;
//    ///   double Shift=7.5;
//    //   double Shift=11.5;
//    double Shift=7.5;
//   return -log(exp(-0.5*(dist-Shift)*(dist-Shift)));
//   ///  return -log(exp(-1.0*(disp(0)-Shift)*(disp(0)-Shift)));

//   //    return 0.1*(disp(0)-6)*(disp(0)-6)+0.8*(disp(0)-6)+2;
//   //  double Shift=3.0;

//   //return -8*dist;
//   }
//   else if (ImpChoice==TUNEDFUNCTION){
//     //    if (PathData.Path.Communicator.MyProc()==1)
//     //      cerr<<"Into tuned action"<<endl;
//     ///    double currShift=8.34;
//     ///    double currShift=9.1;
//     double currShift=6.0;    return -log(exp(-(disp(0)-currShift)*(disp(0)-currShift)));
//     //    return log(10.0)*(-0.24*dist*dist+0.007*dist*dist*dist+0.0006*dist*dist*dist*dist); //-log(1.0/(dist*dist+0.05));
//   }
//    else if (ImpChoice==RETUNEDFUNCTION){
// //   //    //    if (PathData.Path.Communicator.MyProc()==1)
// //   //    //      cerr<<"Into retuned action"<<endl;
// //   //
// //     if (dist>9){
// //   //      //      if (PathData.Path.Communicator.MyProc()==1)
// //   //      //	cerr<<"leaving returned action"<<endl;
// //       return 99999999;
// //     }
// //     else {






//      return log(10.0)*(-0.29*dist*dist+0.0012*dist*dist*dist*dist); //-log(1.0/(dist*dist+0.05));
// //     }
//   }
// //  else if (ImpChoice==RETUNEDFUNCTION){
//     //     //    if (PathData.Path.Communicator.MyProc()==1)
// //     //      cerr<<"Into retuned action"<<endl;

// //     //    if (dist>9){
// //       //      if (PathData.Path.Communicator.MyProc()==1)
// //       //	cerr<<"leaving returned action"<<endl;
// //     //      return 99999999;
// //       //    }
// //     //    else {
// //       //      if (PathData.Path.Communicator.MyProc()==1)
// //       //	cerr<<"now leaving retuned  action"<<endl;
// //    return log(10.0)*(-0.29*dist*dist+0.0012*dist*dist*dist*dist);
// //       //    }
// //  }
//   else {
//     cerr<<"You haven't give a valid choice!"<<endl;
//     assert(1==2);
//   }

//   //  cerr<<"MY shift is "<<shift<<endl;
//   //  return -log(0.01+(dist*dist)*(0.94*exp(-dist*dist)+0.06));

//   /////  return -log(exp(-0.5*(dist-3.5)*(dist-3.5))); //(dist*dist));
//   //    if (PathData.Path.Communicator.MyProc()==1)
//   //      cerr<<"leaving ina  very weird way"<<endl;

// }



double OpenLoopImportanceClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}


string 
OpenLoopImportanceClass::GetName()
{
  return "OpenLoopImportance";
}
