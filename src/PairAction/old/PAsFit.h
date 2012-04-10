#ifndef PA_S_FIT_H
#define PA_S_FIT_H
#include "PAFitBase.h"
#include "../Splines/BicubicSpline.h"
#ifdef MAKE_FIT
#include "../MPI/Communication.h"
#endif

class PAsFitClass : public PairActionFitClass
{
private:
  bool GridsAreMine;
  Array<double,1> Coefs;
  Array<double,1> Pn;
  Array<double,1> sMax;
//   Array<double,1> UsMax;
//   Array<double,1> dUsMax;
#ifdef MAKE_FIT
  CommunicatorClass Comm;
#endif
public:
  Grid *qgrid, *ygrid;
  Array<MultiBicubicSpline,1> Usplines, dUsplines;
  int Order;
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
  PAsFitClass()
  { 
#ifdef MAKE_FIT
    Comm.SetWorld();
#endif
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAsFitClass()
  { if (GridsAreMine){ delete qgrid; delete ygrid; } }
};

#endif
