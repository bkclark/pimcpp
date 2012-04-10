#ifndef AGG_VOLUME_BIAS_MOVE_H
#define AGG_VOLUME_BIAS_MOVE_H

#include "MoleculeMoveBase.h"

class AVBMove: public MolMoveClass{
  double vol, Vin, Vout;
  double Rin, Rout, RCubeDiff;
	double Pbias;
  int numAccepted,numMoves;
  public:
  AVBMove(PathDataClass &myPathData,IOSectionClass outSection); 
  AVBMove(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection)
    {
    }
	void Read(IOSectionClass &moveInput);
	bool AreBound(int slice, int mol1, int mol2);
	bool AreBound(int slice, int mol1, dVec coord);
	dVec GenerateBound(int slice, int ToMol, int FromMol);
	dVec GenerateUnBound(int slice, int ToMol, int FromMol);
	double Sample(int& slice1, int& slice2, Array<int,1>& activeParticles);
/*  void WriteRatio()
    {
    };*/

};

#endif
