#ifndef WORM_STAGE_CLASS_H
#define WORM_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"

class WormStageClass : public LocalStageClass
{
public:
  void WriteRatio();
  ObservableDouble AcceptRatioVar;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  void Read(IOSectionClass &io);

  void Accept();
  void Reject();
  bool PadWorm();
  void MakeRoomForGrowingTail(int changeAmount,int &tailSlice);
  void MakeRoomForShrinkingTail(int changeAmount,int &tailSlice);
  void MakeRoomForGrowingHead(int slicesToGrow,int &headSlice);
  void MakeRoomForShrinkingHead(int slicesToShrink,int &headSlice);
  bool MoveHead;
  bool Grow; 
  int ChangeAmount;
  
  WormStageClass(PathDataClass &pathData, 
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection),
    AcceptRatioVar("AcceptRatio",OutSection,pathData.Path.Communicator) 
  { 
    //do nothing for now
    BisectionLevel = 0;
  }
};

#endif
