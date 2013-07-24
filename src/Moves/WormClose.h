#ifndef WORM_CLOSE_MOVE_H
#define WORM_CLOSE_MOVE_H


#include "../PathDataClass.h"
#include "EmptyStage.h"
#include "LeviFlightStage.h"
#include "MultiStage.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormCloseMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  LeviFlightStageClass LeviFlightStage;
  EmptyStageClass EmptyStage;
  ///We distinguish between the ActiveParticles which are the one's
/// the stage actually  works on and whose action may change and the
/// SavedParticles which are the particles that need to be saved.
  Array<int,1> SavedParticles;
  int SavedSlice1;
  int SavedSlice2;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  void MakeMove();
  void ChooseTimeSlices();
  void SwapPermutations(int newHead,int oldHead);
  void Accept();
  void Reject();
  // Actually attempts the move and accepts or rejects
  //  void MakeMove();
  WormCloseMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    LeviFlightStage(pathData, outSection),
    EmptyStage(pathData,0,outSection)
    
  {
    DumpFreq = 20;
  }
};


#endif
