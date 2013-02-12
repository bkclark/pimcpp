/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#ifndef LOOP_CLASS_H
#define LOOP_CLASS_H

#include "EventClass.h"

class ObservableClass;
class MoveClass;

/// Note:  LoopClass is not allowed to write any output to its out IOsection.
class LoopClass : public EventClass
{
protected:
  int NumSteps;
  bool Equilibrate;
  std::list<MoveClass*> &Moves;
  std::list<ObservableClass*> &Observables;
  std::list<EventClass*> Events;
  MoveClass *FindMove(string name);
  ObservableClass *FindObservable(string name);
 public:
  void DoEvent();
  void Read(IOSectionClass &IO);
  void Read(IOSectionClass &IO, int steps);
  void SetOutfile (IOSectionClass &out)
  { IOSection = out; }
  LoopClass(PathDataClass &pathData, IOSectionClass &out, list<MoveClass*> &moves, list<ObservableClass*> &observables) :
    EventClass(pathData, out), NumSteps(1), Moves(moves), Observables(observables)
  {
    //  do nothing
  }
};



#endif
