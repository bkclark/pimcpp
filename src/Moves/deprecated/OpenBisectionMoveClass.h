#ifndef OPENBISECTIONMOVE_CLASS2_H
#define OPENBISECTIONMOVE_CLASS2_H

//#include "../Common.h"
#include "MoveBase.h"
#include "../PathDataClass.h"
//#include "../ActionClass.h"
#include "BisectionClass.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class OpenBisectionMoveClass : public ParticleMoveClass
{
 private:
  BisectionClass Bisection;
  string Open;
 public:
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass &moveInput);
  /// Function to actually make a bisection move.
  void MakeMove();

  OpenBisectionMoveClass(PathDataClass &myPathData, IOSectionClass iosection) :
    ParticleMoveClass(myPathData, iosection), 
    Bisection(myPathData)
  { 
    /* Do nothing for now. */ 
  }
};


#endif
