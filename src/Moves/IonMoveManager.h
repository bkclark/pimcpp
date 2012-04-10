#ifndef ION_MOVE_MANAGER_CLASS_H
#define ION_MOVE_MANAGER_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "../Observables/ObservableVar.h"
#include "IonMoveStage.h"

class IonMoveClass : public MultiStageClass
{
private:
  int StepNum;
  int NumLevels, LowestLevel;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  StageClass* MoveStage;
  //  ObservableDouble AcceptanceRatioVar;
public:
  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  //void MakeMove();

  IonMoveClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out),StepNum(0)

  { 
    // do nothing for now
  }
};


#endif
