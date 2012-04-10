#ifndef MOLECULE_MOVE_MANAGER_CLASS_H
#define MOLECULE_MOVE_MANAGER_CLASS_H


#include "../PathDataClass.h"
#include "MoveBase.h"
#include "../Observables/ObservableVar.h"
#include "MoleculeMove.h"
#include "AVBMove.h"
#include "MultiStage.h"
#include "MoleculeBias.h"

class MoleculeMoveStageManagerClass : public MultiStageClass
{
private:
  int SpeciesNum;
//  void ChooseTimeSlices();
  //StageClass* MoveStage;
  StageClass* MoveStage;
  //  ObservableDouble AcceptanceRatioVar;
public:
  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  //void MakeMove();

  MoleculeMoveStageManagerClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out)

  { 
    // do nothing for now
  }
};


#endif
