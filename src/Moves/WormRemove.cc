#include "WormRemove.h"

void
WormRemoveMoveClass::Read (IOSectionClass &in)
{
  // Construct action list
  //LeviFlightStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  //  LeviFlightStage.Read(in);
  //  LeviFlightStage.Actions.push_back(&PathData.Actions.Kinetic);
  //  LeviFlightStage.Actions.push_back(&PathData.Actions.Mu);
  
  //  Stages.push_back(&LeviFlightStage);
  
  //  EmptyStage.Actions.push_back(&PathData.Actions.ShortRange);
  EmptyStage.Actions.push_back(&PathData.Actions.Mu);
  Stages.push_back(&EmptyStage);
  ActiveParticles.resize(1);
}


int 
WormRemoveMoveClass::CountWormSize()
{
  int headSlice;
  int headPtcl;
  int tailSlice;
  int tailPtcl;
  PathData.FindHead(headSlice,headPtcl);
  PathData.FindTail(tailSlice,tailPtcl);
  int currSlice=headSlice;
  int currPtcl=headPtcl;
  int wormSize=0;
  int nextSlice;
  int nextPtcl;


  while (currSlice!=tailSlice && currPtcl!=tailPtcl){

    if (currSlice==PathData.Join)
      nextPtcl=PathData.Path.Permutation(nextPtcl);
    else
      nextPtcl=currPtcl;
    nextSlice=(nextSlice+1) % PathData.Path.NumTimeSlices();
    currSlice=nextSlice;
    currPtcl=nextPtcl;
    wormSize++;
  }
  return wormSize;
}



dVec WormRemoveMoveClass::RandomBoxLocation()
{
  dVec myLoc;
  dVec boxSize=PathData.Path.GetBox();
  myLoc[0]=PathData.Path.Random.Local()*(2*boxSize[0])-boxSize[0];
  myLoc[1]=PathData.Path.Random.Local()*(2*boxSize[1])-boxSize[1];
  myLoc[2]=PathData.Path.Random.Local()*(2*boxSize[2])-boxSize[2];
  return myLoc;
}
 


void 
WormRemoveMoveClass::MakeMove()
{
#if 1==1
  cerr<<"Num Into worm remove class"<<endl;
  Slice1=0;
  Slice2=0;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;

  int headSlice;
  int tailSlice;
  int headPtcl;
  int tailPtcl;
  int currPtcl;
  int wormSize=-1;
  int numEmpty;
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  if (PathData.Path.NowOpen){ 
    PathData.WormInfo(headSlice,headPtcl,
		      tailSlice,tailPtcl,
		      numEmpty,wormSize);

  cerr<<"Worm remove size is "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;

    if (wormSize==0){
      PathData.Path.PrintRealSlices();
      cerr<<"Slices: "<<headSlice<<" "<<headPtcl<<" "
	  <<tailSlice<<" "<<tailPtcl<<endl;
      assert(wormSize!=0);
    }
    if (wormSize!=1){
      MultiStageClass::Reject();
      return;
    }
    if (headPtcl!=tailPtcl){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      PathData.MoveJoin(PathData.NumTimeSlices()-1);
    }
    PathData.WormInfo(headSlice,headPtcl,
		      tailSlice,tailPtcl,
		      numEmpty,wormSize);
    assert(headPtcl==tailPtcl);
    for (int slice=headSlice;slice<=tailSlice;slice++){
      PathData.Path.ParticleExist(slice,headPtcl)=0.0;
    }
    ActiveParticles.resize(1);
    ActiveParticles(0)=headPtcl;
    Slice1=headSlice;
    Slice2=tailSlice;
    PathData.Path.NowOpen=false;
  }
  else {
    cerr<<"Num size worm is 0!"<<endl;
    int randShift=
      PathData.Path.Random.LocalInt(PathData.Path.NumTimeSlices()-2);
    PathData.MoveJoin(0);
    PathData.ShiftData(randShift);
    PathData.Join=randShift;
    PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
    int emptyPtcl=PathData.Path.FindEmptyParticle();
    assert(emptyPtcl!=-1);
    double lambda=PathData.Path.ParticleSpecies(emptyPtcl).lambda;
    double tau=PathData.Path.tau;
    double sigma2=(2.0*lambda*tau);
    double sigma=sqrt(sigma2);

    dVec delta;
    PathData.Path.Random.LocalGaussianVec(sigma,delta);
    PathData.Path(1,emptyPtcl)=RandomBoxLocation();
    PathData.Path(2,emptyPtcl)=PathData.Path(1,emptyPtcl)+delta;
    PathData.Path.ParticleExist(1,emptyPtcl)=1.0;
    PathData.Path.ParticleExist(2,emptyPtcl)=1.0;
    ActiveParticles.resize(1);
    ActiveParticles(0)=emptyPtcl;
    PathData.Path.NowOpen=true;
    Slice1=1;
    Slice2=2;
    
  }
    PathData.WormInfo(headSlice,headPtcl,
		      tailSlice,tailPtcl,
		      numEmpty,wormSize);

  cerr<<"Worm remove size as "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;
  if (numEmpty<PathData.Path.NumTimeSlices()){
    MultiStageClass::Reject();
    return;
  }
  MultiStageClass::MakeMove();
    PathData.WormInfo(headSlice,headPtcl,
		      tailSlice,tailPtcl,
		      numEmpty,wormSize);

  cerr<<"Worm remove size2 as "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;

  cerr<<"Num Out of worm remove class"<<endl;    
#endif      
}
