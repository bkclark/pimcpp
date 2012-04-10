/// \cond
//#include "../PathDataClass.h"
#include "BisectionMoveClass.h"
//#include "../Common.h"
#include "../SpeciesClass.h"


void BisectionMoveClass::Read(IOSectionClass &moveInput)
{
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="Bisection");
  assert(moveInput.ReadVar("name",Name));
  Array<int,1> tempActiveSpecies(0);
  assert(moveInput.ReadVar("ActiveSpecies",tempActiveSpecies));
  SetActiveSpecies(tempActiveSpecies);
  int tempNumParticlesToMove;
  assert(moveInput.ReadVar("NumParticlesToMove",tempNumParticlesToMove));
  SetNumParticlesToMove(tempNumParticlesToMove);
  moveInput.ReadVar("StartTimeSlice",StartTimeSlice);
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


void BisectionMoveClass::WriteRatio()
{
  Array<double,1> ratios;
  Bisection.AcceptanceRatios(ratios);
  if (FirstTime) {
    FirstTime = false;
    Array<double,2> writeRatios(1,ratios.size());
    writeRatios(0,Range::all()) = ratios;
    OutSection.WriteVar("AcceptRatio", writeRatios);
    IOVar = OutSection.GetVarPtr("AcceptRatio");
  }
  else {
    IOVar->Append(ratios);
  }
  OutSection.FlushFile();
}




void BisectionMoveClass::MakeMove()
{
  int numSlices=PathData.NumTimeSlices();
  double xi=PathData.Path.Random.Local();
  int slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
  int slice1=slice2-(1<<NumLevels);
  PathData.MoveJoin(slice2);

  //  int EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
  StartTimeSlice=slice1;
  int EndTimeSlice=slice2;

  //  int EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
  
  if (PathData.Path.OpenPaths){
    if (StartTimeSlice<=(int)PathData.Path.OpenLink && (int)PathData.Path.OpenLink<=EndTimeSlice)
      ChooseParticles(); 
    else 
      ChooseParticlesOpen();
  }
  else {
    ChooseParticlesOpen();
  }
  PathData.MoveJoin(EndTimeSlice);
  bool toAccept=Bisection.Bisect(StartTimeSlice,NumLevels,ActiveParticles);
  if (toAccept ==true ){
    PathData.AcceptMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
    NumAccepted++;
  }
  else {
    PathData.RejectMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
  }
  NumMoves++;
  MoveClass::MakeMove();
}
/// \endcond
