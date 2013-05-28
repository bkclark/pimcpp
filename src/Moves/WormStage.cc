
#include "WormStage.h"

void WormStageClass::Read(IOSectionClass &in)
{

  
}

void WormStageClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

void WormStageClass::Accept()
{
  //do nothing for now
  
}

void WormStageClass::Reject()
{
  //do nothing for now

}

void WormStageClass::MakeRoomForGrowingHead(int slicesToGrow,int &headSlice)
{
  int lastSlice=PathData.Path.NumTimeSlices()-1;
  PathData.MoveJoin(0);
  if (1+slicesToGrow<=headSlice){
    PathData.MoveJoin(lastSlice);
    return;
  }
  int shiftNeeded=1+slicesToGrow-headSlice;
  PathData.Path.ShiftData(shiftNeeded);
  PathData.Join=shiftNeeded;
  PathData.MoveJoin(lastSlice);
  headSlice=headSlice+shiftNeeded;
}


///Because this function pushes the head to a place hwere it can
///shrink slicesToShrink without removing the entire particle from
///this slice we must have that slicesToShrink<NumTimeSlices()-1
void WormStageClass::MakeRoomForShrinkingHead(int slicesToShrink,int &headSlice)
{
  int lastSlice=PathData.Path.NumTimeSlices()-1;
  PathData.MoveJoin(lastSlice);
  if (headSlice+slicesToShrink<lastSlice)
    return;
  int shiftNeeded=headSlice+slicesToShrink-lastSlice+1;
  PathData.Path.ShiftData(-shiftNeeded);
  PathData.Join=lastSlice-shiftNeeded;
  PathData.MoveJoin(lastSlice);
  headSlice=headSlice-shiftNeeded;
}


void WormStageClass::MakeRoomForGrowingTail(int slicesToGrow,int &tailSlice)
{
  int lastSlice=PathData.Path.NumTimeSlices()-1;
  PathData.MoveJoin(lastSlice);
  if (lastSlice-tailSlice>slicesToGrow)
    return;
  int shiftNeeded=slicesToGrow+tailSlice-lastSlice+1;
  PathData.Path.ShiftData(-shiftNeeded);
  PathData.Join=lastSlice-shiftNeeded;
  PathData.MoveJoin(lastSlice);
  tailSlice=tailSlice-shiftNeeded;

}


void WormStageClass::MakeRoomForShrinkingTail(int slicesToShrink,int &tailSlice)
{
  //  cerr<<"My number of time slices is "<<PathData.Path.NumTimeSlices()<<endl;
  int lastSlice=PathData.Path.NumTimeSlices()-1;
  if (tailSlice-slicesToShrink>0){
    PathData.MoveJoin(lastSlice);
    return;
  }
  //    PathData.Path.PrintRealSlices();
  int shiftNeeded=slicesToShrink-tailSlice+1;
  //  cerr<<"Shift needed: "<<shiftNeeded<<endl;
  PathData.MoveJoin(0);
  PathData.Path.ShiftData(shiftNeeded);
  PathData.Join=shiftNeeded;
  PathData.MoveJoin(lastSlice);
  tailSlice=tailSlice+shiftNeeded;
  //  cerr<<"posT"<<endl;
}


bool WormStageClass::PadWorm()
{
  int headSlice;
  int headPtcl;
  int tailSlice;
  int tailPtcl;
  PathData.FindHead(headSlice,headPtcl);
  PathData.FindTail(tailSlice,tailPtcl);
  

  int nextPtcl=tailPtcl;
  int nextSlice=tailSlice;
  int numEmpty=0;
  do{
    if (nextSlice==PathData.Join)
      nextPtcl=PathData.Path.Permutation(nextPtcl);
    nextSlice=(nextSlice+1) % PathData.Path.NumTimeSlices();
    //    if (nextSlice!=0) //We should only count iether the beginning or the end slice but not both of them
      numEmpty++;
  }
  while (PathData.Path.ParticleExist(nextSlice,nextPtcl)==0.0);
  numEmpty--;
  cerr<<"Num empty is "<<numEmpty<<endl;
  //  cerr<<"Num full is "<<PathData.Path.NumTimeSlices()*PathData.Path.NumParticles()-numEmpty-50<<endl;
  if (numEmpty<PathData.Path.NumTimeSlices()+1){
    return true;
  }
  else 
    return false;
//   if (numEmpty<PathData.Path.NumParticles()){
//     int oldNumParticles=PathData.Path.NumParticles();
//     int newNumParticles=oldNumParticles+1;
//     PathData.Path.resizeAndPreserve(PathData.Path.NumTimeSlices(),
// 				    newNumParticles);
//     PathData.Path.ParticleExist.resizeAndPreserve(PathData.Path.NumTimeSlices(),
// 						  newNumParticles);
//     PathData.Path.Permutation.resizeAndPreserve(newNumParticles);
//   }
}

///Because this stage does shifting, it must be called before any
///tentative moves to the path have been called
double WormStageClass::Sample(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  //  cerr<<"hmm"<<endl;
  if (!PathData.Path.NowOpen)
    return 1.0;
  //  cerr<<"Sampling now"<<endl;
  //  cerr<<"The join is "<<PathData.Join<<endl;

  PathClass &Path=PathData.Path;
  int headSlice;
  int headPtcl;
  int tailSlice;
  int tailPtcl;
  int numEmpty;
  int wormSize;
  PathData.WormInfo(headSlice, headPtcl,
			 tailSlice, tailPtcl,
			 numEmpty,  wormSize);
  


  SetMode(OLDMODE);
  //  PathData.Path.PrintRealSlices();
  SetMode(NEWMODE);
  //  cerr<<"The join is "<<PathData.Join<<endl;
  //  cerr<<"My head is at "<<headSlice<<" "<<headPtcl<<endl;
  //  cerr<<"My tail is at "<<tailSlice<<" "<<tailPtcl<<endl;

  double lambda=PathData.Path.ParticleSpecies(tailPtcl).lambda;
  double tau=PathData.Path.tau;
  double sigma2=(2.0*lambda*tau);
  double sigma=sqrt(sigma2);
  double logSampleProb=0.0;
  //  grow=false;
  if (!MoveHead && Grow){
    cerr<<"Worm Preparing to grow tail "<<ChangeAmount<<endl;
    MakeRoomForGrowingTail(ChangeAmount,tailSlice);
    PathData.WormInfo(headSlice, headPtcl,
	     tailSlice, tailPtcl,
	     numEmpty, wormSize);

    for (int slice=tailSlice;slice<tailSlice+ChangeAmount;slice++){
      dVec delta;
      Path.Random.LocalGaussianVec(sigma,delta);
      double dist2=dot(delta,delta);
      logSampleProb+=-0.5*dist2/sigma2;
      PathData.Path(slice+1,tailPtcl)=PathData.Path(slice,tailPtcl)+delta;
      PathData.Path.ParticleExist(slice+1,tailPtcl)=1.0;
    }
    slice1=tailSlice;
    slice2=tailSlice+ChangeAmount;
    activeParticles.resize(1);
    activeParticles(0)=tailPtcl;
    //    cerr<<"My log sample prob is "<<logSampleProb<<endl;
    //    if (PadWorm()){
    //      cerr<<"Pad worm has failed";
    //      return -10000000;
    //    }
    //    else
      return exp(-logSampleProb);
  }
  else if (MoveHead && Grow){
    cerr<<"Worm Preparing to grow the head "<<ChangeAmount<<endl;
    MakeRoomForGrowingHead(ChangeAmount,headSlice);
    PathData.WormInfo(headSlice, headPtcl,
	     tailSlice, tailPtcl,
	     numEmpty, wormSize);

    //    cerr<<headSlice<<" "<<ChangeAmount<<endl;
    for (int slice=headSlice;slice>headSlice-ChangeAmount;slice--){
      dVec delta;
      Path.Random.LocalGaussianVec(sigma,delta);
      double dist2=dot(delta,delta);
      logSampleProb+=-0.5*dist2/sigma2;
      PathData.Path(slice-1,headPtcl)=PathData.Path(slice,headPtcl)+delta;
      PathData.Path.ParticleExist(slice-1,headPtcl)=1.0;
    }
    slice2=headSlice;
    slice1=headSlice-ChangeAmount;
    activeParticles.resize(1);
    activeParticles(0)=headPtcl;
    //    cerr<<"My log sample prob is "<<logSampleProb<<endl;
    //    if (PadWorm()){
    //      cerr<<"I'm rejecting this"<<endl;
    //      return -10000000;
    //    }
    //    else
      return exp(-logSampleProb);
  }
  else if (MoveHead && !Grow){
    cerr<<"Worm Preparing to shrink head "<<ChangeAmount<<endl;
    MakeRoomForShrinkingHead(ChangeAmount,headSlice);
    PathData.WormInfo(headSlice, headPtcl,
	     tailSlice, tailPtcl,
	     numEmpty, wormSize);

    int slice=headSlice;
    slice2=headSlice;
    slice1=headSlice;
    //    cerr<<headSlice<<" "<<tailSlice<<" "<<ChangeAmount<<endl;
    while (slice<headSlice+ChangeAmount && (!(headPtcl==tailPtcl && tailSlice==slice))){
      PathData.Path.ParticleExist(slice,headPtcl)=0.0;
      dVec oldDelta=PathData.Path(slice,headPtcl)-PathData.Path(slice+1,headPtcl);
      PathData.Path.PutInBox(oldDelta);
      double dist2=dot(oldDelta,oldDelta);
      logSampleProb+=-0.5*dist2/sigma2;
      slice++;
      slice2++;
    }
    //    cerr<<slice2<<" "<<headSlice+ChangeAmount<<endl;
    activeParticles.resize(1);
    activeParticles(0)=headPtcl;
    //    if (slice2!=headSlice+changeAmount ||
    //	(!(headPtcl==tailPtcl && tailSlice==slice2))){
    //      cerr<<"I SHOULD BE REJECTING"<<endl;
    //      return -10000000;
    //    }
    //    PadWorm();
    return exp(logSampleProb);
  }
  else if (!MoveHead && !Grow){
    //    PathData.Path.PrintRealSlices();
    cerr<<"Worm Preparing to shrink "<<ChangeAmount<<" "<<tailSlice<<" "
	<<tailPtcl<<" "
	<<headSlice<<" "
	<<headPtcl<<endl;
    MakeRoomForShrinkingTail(ChangeAmount,tailSlice);
    //    PathData.Path.PrintRealSlices();
    PathData.WormInfo(headSlice, headPtcl,
	     tailSlice, tailPtcl,
	     numEmpty, wormSize);

    int slice=tailSlice;
    slice1=tailSlice;
    while (slice>tailSlice-ChangeAmount &&(!(headPtcl==tailPtcl && headSlice==slice))){
      PathData.Path.ParticleExist(slice,tailPtcl)=0.0;
      dVec oldDelta=PathData.Path(slice,tailPtcl)-PathData.Path(slice-1,tailPtcl);
      double dist2=dot(oldDelta,oldDelta);
      logSampleProb+=-0.5*dist2/sigma2;
      slice--;
      slice1--;
    }
    //    slice1=tailSlice-ChangeAmount;
    slice2=tailSlice;
    activeParticles.resize(1);
    activeParticles(0)=tailPtcl;
    //    PathData.Path.PrintRealSlices();
    //    if (slice1!=tailSlice-ChangeAmount || 
    //	(!(headPtcl==tailPtcl && headSlice==slice1))){
    //      cerr<<"I SHOULD BE REJECTING"<<endl;
    //      return -10000000;
    //    }
    //    PadWorm();
    return exp(logSampleProb);
  }
}


//   if (moveHead && grow){
//     int changeCount=ChangeAmount;
//     int newHeadLoc=currHeadLoc;
//     int newHeadPtcl=currHeadPtcl;
//     while (changeCount>0){
//       changeCount--;
//       newHeadLoc++;
//       if (newHeadLoc==currTailLoc && newHeadPtcl==currTailPtcl){
// 	//	Path.Path.resizeAndPreserve(
// 	Path(
	
    
    




//   int skip = 1<<(BisectionLevel+1);
//   double levelTau = 0.5*PathData.Path.tau*skip;

//   int numImages = PathData.Actions.NumImages;

//   double logSampleProb=0.0;
//   double logOldSampleProb=0.0;

//   for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
//     int ptcl=activeParticles(ptclIndex);
//     double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
//     for (int slice=slice1;slice<slice2;slice+=skip){
//       SetMode(OLDMODE);
//       dVec rOld=Path(slice,ptcl);
//       dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
//       dVec rbarOld=rOld+ 0.5*rdiffOld;
//       dVec rppOld=Path(slice+(skip>>1),ptcl);
//       dVec DeltaOld=rppOld-rbarOld;
//       Path.PutInBox(DeltaOld);
      
//       SetMode(NEWMODE);
//       dVec r=Path(slice,ptcl);
//       dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
  
//       dVec rbar=r+ 0.5*rdiff;
//       ///Here we've stored the new position in the path
//       dVec  Delta;
//       Path.Random.LocalGaussianVec(sigma,Delta);
//       PathData.Path.PutInBox(Delta);

//       double GaussProd=1.0;
//       double GaussProdOld=1.0;
//       for (int dim=0; dim<NDIM; dim++) {
// 	double GaussSum = 0.0;
// 	double GaussSumOld =0.0;
// 	for (int image=-numImages; image <= numImages; image++) {
// 	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
// 	  double distOld=DeltaOld[dim]+(double)image*Path.GetBox()[dim];
// 	  GaussSum += exp(-0.5*dist*dist/sigma2);
// 	  GaussSumOld += exp(-0.5*distOld*distOld/sigma2);
// 	}
// 	GaussProd *= GaussSum;
// 	GaussProdOld *= GaussSumOld;
//       }
//       logOldSampleProb += prefactorOfSampleProb + log(GaussProdOld);
//       logSampleProb += prefactorOfSampleProb  + log(GaussProd);
      
//       dVec rpp=rbar+Delta;
//       ////REALLY BAD HACK!
//       ///      dVec rpp=rbar+100*Delta;
//       ///Here we've stored the new position in the path
//       Path.SetPos(slice+(skip>>1),ptcl,rpp);
//     }
//   }
// //   SetMode (NEWMODE);
// //   SamplePaths (slice1, slice2, activeParticles);
// //   SetMode(OLDMODE);
// //   double oldSample = LogSampleProb (slice1, slice2, activeParticles);
// //   SetMode(NEWMODE);
// //   double newSample = LogSampleProb (slice1, slice2, activeParticles);
//   //cerr << "oldSample = " << oldSample << endl;
//   //cerr << "oldSampleProb = " << logOldSampleProb << endl;
//   //cerr << "newSample = " << newSample << endl;
//   //cerr << "newSampleProb = " << logSampleProb << endl;

//   return exp(-logSampleProb+logOldSampleProb);
//   //return (exp (-newSample + oldSample));
// }

