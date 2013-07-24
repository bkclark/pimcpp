#ifndef DIAGONAL_ACTION_CLASS
#define DIAGONAL_ACTION_CLASS
#include "ActionBase.h"
#include "../PairAction/PAFit.h"



/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class DiagonalActionClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;

  int TotalTime;
  /// These are the coefficients used for the low-variance estimator
  /// for the gradient
  Array<double,1> ck;
  int NumBasisFuncs, m;
  double Router;
  bool UseLowVariance;
  void Setup_ck();
  inline double g(double r);
public:

  
  bool HaveSamplingTable;
  void Read (IOSectionClass &in);
  /*   double dUdR(int slice,int ptcl1, int ptcl2, int level); */
  /*   double  dUdR_movers(int slice,int ptcl1, int ptcl2, int level); */
  /*   double d2UdR2(int slice,int ptcl1, int ptcl2, int level); */
  /*   double d2UdR2_movers(int slice,int ptcl1, int ptcl2, int level); */
  

  
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  /*   double SingleActionForcedPairAction (int slice1, int slice2, */
  /* 				       PairActionFitClass &PA); */
  /*   double  d_dBetaForcedPairAction (int slice1, int slice2, */
  //    PairActionFitClass &pA);
  
  double d_dBeta (int slice1, int slice2, int level);
  /*   void GradAction (int slice1, int slice2, const Array<int,1> &ptcls, */
  /* 		   int level, Array<dVec,1> &gradVec); */
  string GetName();
  DiagonalActionClass (PathDataClass &pathData,
		       Array<PairActionFitClass*, 2> &pairMatrix);
};




#endif
