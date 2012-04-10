#ifndef ION_MOVE_STAGE_CLASS_H
#define ION_MOVE_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"

class IonStageClass : public LocalStageClass
{
  Array<int,1> PtclRoster;
public:
  IonStageClass(PathDataClass &pathData, int level, IOSectionClass outSection, string SpeciesName); 
  int specID;
  double dRMag;
  void Set(double setdRMag);
  void WriteRatio();
  void Read(IOSectionClass& in);
	bool Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange);
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
//  void Accept();
//  void Reject();
  //int TotalLevels;
};

#endif
