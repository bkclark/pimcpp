#ifndef WORM_REMOVE_MOVE_H
#define WORM_REMOVE_MOVE_H


#include "../PathDataClass.h"
#include "EmptyStage.h"
#include "MultiStage.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormRemoveMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  EmptyStageClass EmptyStage;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  void MakeMove();
  dVec RandomBoxLocation();
  int CountWormSize();
  // Actually attempts the move and accepts or rejects
  //  void MakeMove();
  WormRemoveMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    EmptyStage(pathData,0,outSection)

  {
    DumpFreq = 20;
  }
};


#endif
