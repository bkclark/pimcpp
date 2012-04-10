#include "PathClass.h"


void PathClass::MoveJoinParticleExist(int oldJoin, int newJoin)
{
  bool swappedAlready=false;
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
	ParticleExist[OLDMODE](timeSlice,ptcl)=ParticleExist[NEWMODE](timeSlice,Permutation(ptcl));
      //Now that we've copied the data from B into A, we need to copy the 
      //information into B
      for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){ 
	for (int ptcl=0;ptcl<NumParticles();ptcl++){
	  ParticleExist[NEWMODE](timeSlice,ptcl)=ParticleExist[OLDMODE](timeSlice,ptcl);
	  
	}
      }
    }
  }
  else if (oldJoin>newJoin){
    //  else if (oldJoin>=newJoin){//CHANGED!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	ParticleExist[OLDMODE](timeSlice,Permutation(ptcl))=ParticleExist[NEWMODE](timeSlice,ptcl);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
	ParticleExist[NEWMODE](timeSlice,ptcl)=ParticleExist[OLDMODE](timeSlice,ptcl);
    }
  }
}



void PathClass::ShiftParticleExist(int slicesToShift)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numPtcls=NumParticles()+OpenPaths;
  int numSlices=NumTimeSlices();
  //  cerr<<"SLICES TO SHIFT: "<<slicesToShift<<" "<<numSlices<<endl;
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
	ParticleExist[NEWMODE](slice,ptcl) = ParticleExist[NEWMODE](slice-slicesToShift,ptcl);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
	ParticleExist[NEWMODE](slice,ptcl) = ParticleExist[NEWMODE](slice-slicesToShift,ptcl);
  }
  

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numPtcls;
  Array<double,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0){
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
	///If shifting forward, don't send the last time slice (so always)
	///send slice-1
	sendBuffer(buffIndex)=ParticleExist[OLDMODE](slice-1,ptcl);
	buffIndex++;
      }
    }
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
	///If shifting backward, don't send the first time slice (so always)
	///send slice+1
	sendBuffer(buffIndex)=ParticleExist[OLDMODE](slice+1,ptcl);
	buffIndex++;
      }
    }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
  if (slicesToShift>0)
    startSlice=0;
  else 
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
    for (int ptcl=0;ptcl<numPtcls;ptcl++){
      ParticleExist[NEWMODE](slice,ptcl)=receiveBuffer(buffIndex);
      buffIndex++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int ptcl=0; ptcl<numPtcls; ptcl++)
      ParticleExist[OLDMODE](slice,ptcl) = ParticleExist[NEWMODE](slice,ptcl);

  //  cerr<<"Out of shiftexist"<<endl;
}




  ///////////////////////////
  //// Worm Moves        ///
  /////////////////////////
int PathClass::FindEmptyParticle()
{
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    if (IsEmpty(ptcl))
      return ptcl;
  }
  return -1;
}

void PathClass::FindHead(int &headSlice,int &headPtcl)
{
  
  //do nothing for now
}
void PathClass::FindTail(int &tailSlice,int &tailPtcl)
{
  //do nothing for now
}


bool PathClass::IsEmpty(int ptcl)
{
  for (int slice=0;slice<NumTimeSlices();slice++)
    if (ParticleExist(slice,ptcl)!=0)
      return false;
  return true;

}


bool PathClass::IsFull(int ptcl)
{
  for (int slice=0;slice<NumTimeSlices();slice++)
    if (ParticleExist(slice,ptcl)!=1)
      return false;
  return true;

}


//void MaxAllowedParticles()
//{
//  return 
//
//}

//Initializes the array that specifies what the real slices are.
//Must be called after the number of particles is set
void PathClass::InitRealSlices()
{

  //  return;
  SetMode(NEWMODE);
//   //  RealSlices.resize(MaxAllowedParticles());
//   RealSlices.resize(NumParticles());
//   for (int ptcl=0;ptcl<RealSlices.size();ptcl++){
//     RealSlices(ptcl).first=0;
//     RealSlices(ptcl).last=NumTimeSlices()-1;
//   }
//   RealSlices.AcceptCopy();
  for (int slice=0;slice<NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<NumParticles();ptcl++)
      ParticleExist(slice,ptcl)=1.0;

  if (!WormOn){ //if the worm is not on then all the slices are
		//actually real.  Strictly speaking, this  should
		//really not matter because it isn't being used.
    ParticleExist.AcceptCopy();
    return;
  }
  cerr<<"Worm Initialization"<<endl;
  NowOpen=true;
  NowOpen.AcceptCopy();
  ///If you are here then you want to be initializing the worm.
  ///Specifies how many empty "extra" particles should you start with.

  //  assert(in.ReadVar("StartEmpty",startEmpty));
  int startEmpty=2;
  for (int ptcl=NumParticles()-startEmpty;ptcl<NumParticles();ptcl++)
    for (int slice=0;slice<NumTimeSlices();slice++){
      ParticleExist(slice,ptcl)=0.0;
    }
  ParticleExist(NumTimeSlices()-1,NumParticles()-startEmpty-1)=0.0;
  ParticleExist(0,NumParticles()-startEmpty-1)=0.0;
  ParticleExist.AcceptCopy();
  
  for (int ptcl=NumParticles()-startEmpty-1;ptcl<NumParticles()-1;ptcl++){
    Permutation(ptcl)=ptcl+1;
  }
  Permutation(NumParticles()-1)=NumParticles()-startEmpty-1;
  for (int ptcl=0;ptcl<NumParticles();ptcl++)
    cerr<<"Worm Perm"<<Permutation(ptcl)<<endl;
  Permutation.AcceptCopy();
  Path.AcceptCopy();
}

void PathClass::PrintRealSlices()
{
  for (int ptcl=0;ptcl<NumParticles();ptcl++)
    for (int slice=0;slice<NumTimeSlices();slice++){
      if (ParticleExist(slice,ptcl)==1.0)
	cerr<<ptcl<<" "<<slice<<" "<<Path(slice,ptcl)<<endl;
      else 
	cerr<<ptcl<<" "<<slice<<" "<<-1<<endl;
      //      if (RealSlices(ptcl).first<=slice && slice<=RealSlices(ptcl).last)
      //	cerr<<ptcl<<" "<<slice<<" "<<Path(slice,ptcl)<<end;l
      //      else 
      //	cerr<<ptcl<<" "<<slice<<" "<<-1<<endl;
    }

}

//Must be at least three particles allocated to test the real slice shifting code.
void PathClass::TestRealSlices()
{
//   for (int ptcl=0;ptcl<Permutation.size();ptcl++)
//     cerr<<Permutation(ptcl)<<" ";
//   cerr<<endl;
//   InitRealSlices();
  assert(NumParticles()>=3);
  for (int slice=0;slice<NumTimeSlices()-1;slice++){
    Path(slice,0)=slice;
    Path(slice,1)=slice+100;
    Path(slice,2)=slice+200;
  }
  Path(NumTimeSlices()-1,0)=0;
  Path(NumTimeSlices()-1,1)=100;
  Path(NumTimeSlices()-1,2)=200;
  Path.AcceptCopy();
  ParticleExist(0,0)=0;
  ParticleExist(1,0)=0;
  ParticleExist(2,0)=0;
  ParticleExist(NumTimeSlices()-1,0)=0;

//   RealSlices(0).first=3;
//   RealSlices(0).last=8;
//   RealSlices(2).first=-1;
//   RealSlices(2).last=-1;
//   RealSlices(3).first=-1;
//   RealSlices(3).last=-1;
  Permutation(3)=0;  
  Permutation(0)=2;
  Permutation(2)=3;
  Permutation.AcceptCopy();
  //  cerr<<"First printing"<<endl;
  //  PrintRealSlices();
  int shiftAmount=3;
  MoveJoin(NumTimeSlices()-1,0);

  //  ShiftRealSlices(shiftAmount);
  ShiftData(shiftAmount);

  //  cerr<<"Second printing"<<endl;
  //  PrintRealSlices();
}

void PathClass::ShiftRealSlices(int numToShift)
{
  ShiftParticleExist(numToShift);
//   int lastSlice=NumTimeSlices()-1;
//   if (numToShift>0){
//     for (int ptcl=0;ptcl<NumParticles();ptcl++){
//       int emptyA,fullB,emptyC;
//       if (RealSlices(ptcl).first==-1){
// 	emptyA=numToShift;
// 	fullB=0;
// 	emptyC=0;
//       }
//       else{
// 	emptyA=lastSlice-RealSlices(ptcl).last;
// 	fullB=RealSlices(ptcl).last-RealSlices(ptcl).first;
// 	emptyC=RealSlices(ptcl).first-0;
// 	emptyA = (emptyA<numToShift) ? emptyA : numToShift;
// 	fullB = (emptyA+fullB<numToShift) ? fullB : numToShift-emptyA;
// 	emptyC = (emptyC+fullB+emptyA<numToShift) ? emptyC : numToShift-fullB-emptyA;
//       }
//       cerr<<"I'm sending "<<emptyA<<" "<<fullB<<" "<<emptyC<<endl;
//       int myProc=Communicator.MyProc();
//       int numProcs=Communicator.NumProcs();
//       int recvEmptyA;
//       int recvFullB;
//       int recvEmptyC;
//       Communicator.SendReceive (myProc,emptyA,
// 				myProc+1 % numProcs,recvEmptyA);

//       Communicator.SendReceive (myProc,fullB,
// 				myProc+1 % numProcs,recvFullB);

//       Communicator.SendReceive (myProc,emptyC,
// 				myProc+1 % numProcs,recvEmptyC);
//       assert(recvEmptyA+recvFullB+recvEmptyC==numToShift);

//       ///If you are pushing full particles into this particle then you need
//       ///to be totally empty or totally full
//       //      if (!IsEmpty(ptcl) && !IsFull(ptcl) && recvFullB>0){
//       //	cerr<<"If you are shifting data into a particle it must be empty of full"<<endl;
//       //	abort();
//       //      }
//       if (IsEmpty(ptcl) && recvFullB==0){
// 	int dummy=5;
//       }
//       else if (IsEmpty(ptcl)){
// 	RealSlices(ptcl).first=recvEmptyC;
// 	RealSlices(ptcl).last=recvEmptyC+recvFullB;
//       }
//       else if (IsFull(ptcl)) {
// 	RealSlices(ptcl).first=recvEmptyA+recvEmptyC;
// 	RealSlices(ptcl).last=lastSlice;
// 	if (RealSlices(ptcl).last<RealSlices(ptcl).first){ //you are now empty
// 	  RealSlices(ptcl).first=-1;
// 	  RealSlices(ptcl).last=-1;
// 	}
//       }
//       else {
// 	RealSlices(ptcl).first=RealSlices(ptcl).first+recvEmptyA+recvEmptyC;
// 	RealSlices(ptcl).last=RealSlices(ptcl).last+recvEmptyA+recvFullB+recvEmptyC;
// 	if (RealSlices(ptcl).last>lastSlice)
// 	  RealSlices(ptcl).last=lastSlice;
// 	if (RealSlices(ptcl).last<RealSlices(ptcl).first){ //you are now empty
// 	  RealSlices(ptcl).first=-1;
// 	  RealSlices(ptcl).last=-1;
// 	}
//       }
//     }
//   }   
//   else if (numToShift<0){
//     for (int ptcl=0;ptcl<NumParticles();ptcl++){
//       int emptyA=RealSlices(ptcl).first-0; 
//       int fullB=RealSlices(ptcl).last-RealSlices(ptcl).first;
//       int emptyC=lastSlice-RealSlices(ptcl).last;
//     }
//   }

}



///Not worm relevant yet///
///Multistep moves/////

void PathClass::InitializeStateSaving()
{
  ///Must be called after the number of particles and number of timeslices have been initialized
  SavedPath.resize(NumTimeSlices(),NumParticles());
  SavedPermutation.resize(NumParticles());
  
}

//Saves the current state of newmode into an array. This is so you
//could restore it after doing a multistep move or some form of
//correlated sampling
void PathClass::SaveState()
{
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
      for (int slice=0;slice<NumTimeSlices();slice++){
      SavedPath(slice,ptcl)=Path(slice,ptcl);
    }
      SavedPermutation(ptcl)=Permutation(ptcl);
  }
  SavedJoin=Join;

}


//Restores the saved state into the array in whatever mode you are
//currently in 
void PathClass::RestoreState()
{
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    for (int slice=0;slice<NumTimeSlices();slice++){
      Path(slice,ptcl)=SavedPath(slice,ptcl);
    }
    Permutation(ptcl)=SavedPermutation(ptcl);
  }
  Join=SavedJoin;
  
}


