#ifndef MOLECULE_MOVE_BASE_H
#define MOLECULE_MOVE_BASE_H

//#include "MoveBase.h"
#include "StageClass.h"

// Identifies whether moves should be attempted globally
// or one particle at a time (SINGLE)
enum MoveMode{GLOBAL, SINGLE, SEQUENTIAL, MULTIPLE};

//class MolMoveClass: public ParticleMoveClass{
class MolMoveClass: public LocalStageClass{
	int numMol, molIndex;
  int numActionsToRead, startIndex;
  public:
	MoveMode mode;
	int numAccepted, numMoves;
  double rc;
	// located in PathClass now
  // This will store an array of active particles for each molecule in the sim;
  // should eliminate the need to repeatedly asemble these arrays at each move
  //Array<Array<int,1>,1> MolMembers;
	Array<int, 1> MoveList;
  MolMoveClass(PathDataClass&, IO::IOSectionClass);
  //MolMoveClass(PathDataClass&, IO::IOSectionClass, int numToRead, int start);
  dVec GetCOM(int slice, int mol);
  dVec TranslateMol(int slice, Array<int,1>& activePtcls, double epsilon);
  dVec TranslateMol(int slice, Array<int,1>& activePtcls, dVec translate);
  dVec TranslateMolAll(Array<int,1>& activePtcls, double epsilon);
	void TranslatePtcl(int slice, int ptcl, double Sigma);
	void MoveDimerSeparation(int slice, Array<int,1> mol1, Array<int,1> mol2, double Sigma);
  void RotateMol(int slice, Array<int,1>& activePtcls, dVec& axis, double theta);
  void RotateMol(int slice, Array<int,1>& activePtcls, double theta);
  void RotateMolXYZ(int slice, Array<int,1>& activePtcls, double theta);
  void RotateMolXYZ(Array<int,1>& Slices, Array<int,1>& activePtcls, double theta);
  void RotateMolXYZAll(Array<int,1>& activePtcls, double theta);
  void StressAngle(int slice, int ptcl, dVec axis, double theta);
  void StressBond(int slice, int ptcl, int mol, double s);
  void Read (IOSectionClass &in);
	void Accept();
	void Advance();
  void LoadActions(list<ActionBaseClass*> actions);
};

dVec ArbitraryRotate(dVec axis,dVec coord, double phi);
dVec RotateXYZ(dVec coord,int u1,int u2,int u3,double theta);

#endif
