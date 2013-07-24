//#include "../PathDataClass.h"
#include "OpenBisectionMoveClass.h"
//#include "../Common.h"
//#include "../SpeciesClass.h"


void OpenBisectionMoveClass::Read(IOSectionClass &moveInput)
{
  string typeCheck;
  string open;
  assert(moveInput.ReadVar("open",Open));
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="OpenBisection");
  assert(moveInput.ReadVar("name",Name));
  Array<int,1> tempActiveSpecies(0);
  assert(moveInput.ReadVar("ActiveSpecies",tempActiveSpecies));
  SetActiveSpecies(tempActiveSpecies);
  int tempNumParticlesToMove;
  assert(moveInput.ReadVar("NumParticlesToMove",tempNumParticlesToMove));
  SetNumParticlesToMove(tempNumParticlesToMove);
  string tempNumLevels="";
  assert(moveInput.ReadVar("NumLevels",tempNumLevels));
  if (tempNumLevels=="Max"){
    NumLevels=PathData.Action.MaxLevels;
  }
  else {
    cerr<<"Don't know how to different number of levels from max yet"<<endl;
    assert(1==2);
  }

}

///HACK! HACK! HACK!
double sign(double num)
{
  if (num<0)
    return -1.0;
  else return 1.0;
}
    

void OpenBisectionMoveClass::MakeMove()
{

  int StartTimeSlice;
  int EndTimeSlice;
  if (Open=="head"){
    StartTimeSlice=(int)PathData.Path.OpenLink;
    EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
    while (EndTimeSlice>=PathData.Path.NumTimeSlices() || StartTimeSlice==0){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      StartTimeSlice=(int)PathData.Path.OpenLink;
      EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
    }
  }
  else if (Open=="tail"){
    //    cerr<<"Trying to move the tail"<<PathData.Path(PathData.Path.OpenLink,PathData.NumParticles())<<endl;
    EndTimeSlice=(int)PathData.Path.OpenLink;
    StartTimeSlice=EndTimeSlice-(1<<NumLevels);
    while (StartTimeSlice<0 || EndTimeSlice==PathData.Path.NumTimeSlices()-1){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      EndTimeSlice=(int)PathData.Path.OpenLink;
      StartTimeSlice=EndTimeSlice-(1<<NumLevels);
    }



  }
  else {
    cerr<<"ERROR! We don't know whether to do head or tail?"<<endl;
    assert(1==2);
  }

  PathData.MoveJoin(EndTimeSlice);

  ChooseParticles(); 
  ActiveParticles(0)=(int)PathData.Path.OpenPtcl;
  int changePtcl;
  if (Open=="head"){
    //    cerr<<"Setting to head"<<endl;
    changePtcl=(int)PathData.Path.OpenPtcl;
  }
  else if (Open=="tail"){
    //    cerr<<"Setting to Tail"<<endl;
    changePtcl=PathData.Path.NumParticles();
  }
  else {
    cerr<<"I'm neither seeing a head or tail here! Warning! Warning!"<<endl;
  }
  dVec oldPos=PathData.Path(PathData.Path.OpenLink,changePtcl);
  dVec newPos;///was /10 instead of /40 for the free particles
  newPos(0)= oldPos[0]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(0)/40.0)*
    sign(PathData.Path.Random.Local()-0.5);
  newPos(1)= oldPos[1]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(1)/40.0)*
    sign(PathData.Path.Random.Local()-0.5);
  newPos(2)=  oldPos[2]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(2)/40.0)*
    sign(PathData.Path.Random.Local()-0.5);
  if (fabs(newPos(0))>500 || fabs(newPos(1))>500 || fabs(newPos(2))>500){
    cerr<<"ERROR! ERROR! ERROR!"<<newPos(0)<<" "<<newPos(1)<<" "<<newPos(2)<<endl;
  } 
  PathData.Path.SetPos(PathData.Path.OpenLink,changePtcl,newPos);
  //  cerr<<"Current state is "<<changePtcl<<" "<<PathData.NumParticles()<<" "<<newPos<<" "<<PathData.Path(PathData.Path.OpenLink,PathData.NumParticles())<<endl;

  bool toAccept=Bisection.Bisect(StartTimeSlice,NumLevels,ActiveParticles);

  if (toAccept ==true ){
    //    cerr<<"I've accepted the move "<<PathData.Path(PathData.Path.OpenLink,PathData.NumParticles())<<endl;
    PathData.AcceptMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
    NumAccepted++;
  }
  else {
    //    cerr<<"I've rejected the move "<<PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())<<endl;
    PathData.RejectMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
  }
  NumMoves++;
		     
    
}
