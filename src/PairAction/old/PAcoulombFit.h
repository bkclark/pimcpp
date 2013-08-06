#ifndef PA_COULOMB_FIT_H
#define PA_COULOMB_FIT_H

#include "PAFitBase.h"

class PAcoulombFitClass : public PairActionFitClass
{
private:
  bool GridIsMine;
  Array<double,1> Ucoefs;
public:
  int Order;
  Grid *qgrid;
  Array<MultiCubicSpline,1> Ujshort, Ujlong;
  Array<MultiCubicSpline,1> dUjshort, dUjlong;
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
#endif
  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  PAcoulombFitClass()
  { 
    GridIsMine = false; 
    NumBetas=0;
  }
  ~PAcoulombFitClass()
  { if (GridIsMine) delete qgrid;  }
};

#endif
