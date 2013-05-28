#include "WormClose.h"

void
WormCloseMoveClass::Read (IOSectionClass &in)
{
  // Construct action list
  //  LeviFlightStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  LeviFlightStage.Read(in);
  LeviFlightStage.Actions.push_back(&PathData.Actions.Kinetic);
  LeviFlightStage.Actions.push_back(&PathData.Actions.Mu);
  EmptyStage.Actions.push_back(&PathData.Actions.Mu);
  Stages.push_back(&LeviFlightStage);
  Stages.push_back(&EmptyStage);
  
}


    
 


//Assumes that the join is at the last slice
void
WormCloseMoveClass::SwapPermutations(int newHead,int oldHead)
{
  int temp=PathData.Path.Permutation(newHead);
  PathData.Path.Permutation(newHead)=PathData.Path.Permutation(oldHead);
  PathData.Path.Permutation(oldHead)=temp;
}




void WormCloseMoveClass::Accept()
{
  cerr<<"size Accepting"<<endl;
  PathData.AcceptMove(SavedSlice1,SavedSlice2,SavedParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();
       stageIter!=Stages.end();stageIter++){
    (*stageIter)->Accept();
  }  
  NumAccepted++;
}

void WormCloseMoveClass::Reject()
{
  cerr<<"Rejecting"<<endl;
  PathData.RejectMove(SavedSlice1,SavedSlice2,SavedParticles);
  for (list<StageClass*>::iterator 
	 stageIter=Stages.begin();stageIter!=Stages.end();stageIter++){
    (*stageIter)->Reject();
  }
}



void 
WormCloseMoveClass::MakeMove()
{
#if 1==1

  cerr<<"Starting WormCloseClass"<<endl;
  Slice1=0;
  Slice2=0;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;


  bool toAccept=false;
  PathClass &Path=PathData.Path;
  int wormSize,numEmpty;
  int tailSlice,tailPtcl;
  int headSlice,headPtcl;
   PathData.WormInfo(headSlice, headPtcl,
 		    tailSlice, tailPtcl,
 		    numEmpty, wormSize);
  
//   cerr<<"Worm close size is "<<wormSize<<" "<<numEmpty<<" "
//       <<headSlice<<" "<<headPtcl<<" "
//       <<tailSlice<<" "<<tailPtcl<<" "<<endl;

  if (PathData.Path.NowOpen){
    cerr<<"Now open"<<(double)wormSize/PathData.Path.NumTimeSlices()<<endl;
    PathData.MoveTailToSliceZero() ;
    PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
    PathData.WormInfo(headSlice, headPtcl,
		      tailSlice, tailPtcl,
		      numEmpty, wormSize); 
    
    ///asserts to make sure things are ok
    for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
      assert(Path.ParticleExist(slice,tailPtcl)==0);
    }
    /////
    
    if (headSlice>6){
      MultiStageClass::Reject();
      return;
    }

    for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
      Path(slice,tailPtcl)=Path(slice,headPtcl);
      Path.ParticleExist(slice,tailPtcl)=1; 
      Path.ParticleExist(slice,headPtcl)=0;
    }
    
    SwapPermutations(headPtcl,tailPtcl);
    SavedParticles.resize(2);
    SavedParticles(0)=headPtcl;
    SavedParticles(1)=tailPtcl;
    SavedSlice1=0;
    SavedSlice2=PathData.Path.NumTimeSlices()-1;
    ActiveParticles.resize(1);
    ActiveParticles(0)=tailPtcl;
    Slice1=0;
    Slice2=headSlice;
    PathData.Path.NowOpen=false;

    int numFull=0;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
      if (PathData.Path.IsFull(ptcl)){
	numFull++;
      }
    cerr<<"numFull: "<<numFull<<" "<<headSlice<<" "<<PathData.Path.NumTimeSlices()<<endl;
    double prevActionChange=-log((double)numFull)-log((double)PathData.Path.NumTimeSlices()-1)-log(5.0);
    toAccept=LeviFlightStage.Attempt(Slice1, Slice2 ,
				     ActiveParticles,
				     prevActionChange);
    cerr<<"PERM: "
	<<PathData.Path.Permutation(0)<<" "
	<<PathData.Path.Permutation(1)<<" "
	<<PathData.Path.Permutation(2)<<" "<<endl;
  }
  else { ///the worm is not open!

    PathData.Path.NowOpen=true;
    PathData.MoveJoin(0);
    int randShift=PathData.Path.Random.LocalInt(PathData.Path.NumTimeSlices()-2);
    //BAd --->    int randShift=0;
    PathData.ShiftData(randShift);
    PathData.Join=randShift;
    PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
    int emptyPtcl=PathData.Path.FindEmptyParticle();
    assert(emptyPtcl!=-1);

    int numFull=0;
    bool existsFullPtcl=false;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
      if (PathData.Path.IsFull(ptcl)){
	existsFullPtcl=true;
	numFull++;
      }
    cerr<<"Now closed"<<numFull<<endl;
    if (!existsFullPtcl){
      MultiStageClass::Reject();
      return;
    }
    int fullPtcl=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
    PathData.Path.PrintRealSlices();
    while (!PathData.Path.IsFull(fullPtcl))
      fullPtcl=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
    SwapPermutations(emptyPtcl,fullPtcl);
    cerr<<"out of while loop"<<endl;

    SavedParticles.resize(2);
    SavedParticles(0)=emptyPtcl;
    SavedParticles(1)=fullPtcl;
    ActiveParticles.resize(1);
    ActiveParticles(0)=emptyPtcl;
    for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
      PathData.Path(slice,emptyPtcl)=PathData.Path(slice,fullPtcl);
      PathData.Path.ParticleExist(slice,emptyPtcl)=1;
      PathData.Path.ParticleExist(slice,fullPtcl)=0;
    }
    //    int chooseAmountToOpen=PathData.Path.Random.LocalInt(PathData.Path.NumTimeSlices()-1);
    int chooseAmountToOpen=PathData.Path.Random.LocalInt(5)+1;
    SavedSlice1=0;
    SavedSlice2=PathData.Path.NumTimeSlices()-1;
    Slice1=0;
    Slice2=chooseAmountToOpen+1;
    for (int slice=1;slice<=chooseAmountToOpen;slice++)
      PathData.Path.ParticleExist(slice,emptyPtcl)=0;
    
    dVec disp=PathData.Path(chooseAmountToOpen+1,emptyPtcl)-PathData.Path(0,fullPtcl);
    PathData.Path.PutInBox(disp);
    double dist2 = dot(disp,disp);
    int speciesNum=PathData.Path.ParticleSpeciesNum(emptyPtcl);
    
    
    cerr<<"numFull: "<<numFull<<" "<<headSlice<<" "<<PathData.Path.NumTimeSlices()<<endl;    
    double prevActionChange=dist2/(4.0*PathData.Path.Species(speciesNum).lambda*(chooseAmountToOpen+1)*PathData.Path.tau)+log((double)numFull)+log((double)PathData.Path.NumTimeSlices()-1)+log(5.0);
    toAccept=EmptyStage.Attempt(Slice1, Slice2, 
				     ActiveParticles,
				     prevActionChange);
    
  }
  
    if (toAccept){
    cerr<<"PERMO: "
	<<PathData.Path.Permutation(0)<<" "
	<<PathData.Path.Permutation(1)<<" "
	<<PathData.Path.Permutation(2)<<" "<<endl;

      Accept();
    }
    else {
      Reject();
    }
    PathData.Path.PrintRealSlices();
#endif

}

