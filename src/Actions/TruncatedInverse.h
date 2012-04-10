#ifndef TRUNCATED_INVERSE__CLASS_H
#define TRUNCATED_INVERSE__CLASS_H

#include "../MirroredClass.h"
#include "NodalActionClass.h"
#include <Common/Splines/CubicSpline.h>
#include <vector>
struct doubleint
{
  double distance;
  int ptcl;
};


/// VartionalPIClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class TruncatedInverseClass : public NodalActionClass
{
private:
  PathClass &Path;
  double CutoffAverage;
  void calc_u();
  Array<double,1> u;
  Array<double,1> newCol;
  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<double,2> DetMatrix;
  Array<double,2> SmallDetMatrix;
  Array<double,2> SmallDetMatrixOld;
  void BuildOldMatrix(vector<doubleint> determinantPtcl,int size,
		       Array<double,2> &SmallDetMatrix);
  void BuildNewMatrix(vector<doubleint> determinantPtcl,int size,
		      Array<double,2> &SmallDetMatrix);
  double GetDistance(dVec sortFrom,int ptcl2);

  Mirrored1DClass<double> DeterminantList;
  Mirrored2DClass<double> OtherInfo;

public:
  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  string GetName();
  void BuildDeterminantMatrix();
  void BuildSmallDeterminantMatrix();
  double  MinDistance(dVec oldDvec, dVec newDvec, int ptcl);
  void CheckDeterminantMatrix();
  void Read (IOSectionClass &in);
  bool IsGroundState();
  bool IsPositive(int x);
  NodeType Type();
  int NumTimes;
  void AcceptCopy(int slice1, int slice2);
  void RejectCopy(int slice1, int slice2);
  void WriteInfo (IOSectionClass &out);
  int ChangedColumn;
  dVec olddvec;
  dVec newdvec;
  
  TruncatedInverseClass (PathDataClass &pathData);
};

#endif

