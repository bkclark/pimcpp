#include "IonMoveManager.h"

void IonMoveClass::Read(IOSectionClass &in)
{
  string speciesName;
  assert (in.ReadVar ("Species", speciesName));
  MoveStage = new IonStageClass(PathData, NumLevels, IOSection, speciesName);
  MoveStage->Read (in);
  Stages.push_back (MoveStage);
}
