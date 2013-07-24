/// \cond deprecated
#ifndef BISECTIONMOVE_CLASS2_H
#define BISECTIONMOVE_CLASS2_H

//#include "../Common.h"
#include "MoveBase.h"
#include "../PathDataClass.h"
// #include "../ActionClass.h"
#include "BisectionClass.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionMoveClass : public ParticleMoveClass
{
 private:
  BisectionClass Bisection;
 public:
  ///This is the slice in which the bisection move starts.  It ends up
  ///going to StartTimeSlice+2^NumLevels
  int StartTimeSlice; 
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass &moveInput);
  /// Function to actually make a bisection move.
  void MakeMove();
  void WriteRatio();

  BisectionMoveClass(PathDataClass &myPathData, IOSectionClass outSection) : 
    ParticleMoveClass(myPathData, outSection), 
    Bisection(myPathData)
  { 
    ///Defaults to the 0'th time slice but shouldn't matter because it
    ///chooses a different time slice each time.
    StartTimeSlice=0;
  }
};


#endif
/// \endcond
