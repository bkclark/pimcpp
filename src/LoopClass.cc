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

#include "PathDataClass.h"
#include "LoopClass.h"
#include "WriteData.h"
#include "SwitchClass.h"
#include "Moves/MoveBase.h"
#include "Observables/ObservableBase.h"
#include <sstream>
#include <sys/time.h>

void
LoopClass::DoEvent()
{
  PathData.Path.Equilibrate = Equilibrate;
  std::list<EventClass*>::iterator iter;
  for (int step=0; step<NumSteps; step++){
    for (iter=Events.begin(); iter!=Events.end(); iter++) {
     stringstream output;
     output<<(*iter)->Name<<"("<<step<<"/"<<NumSteps<<")("<<PathData.Path.CloneStr<<") "<<endl;
     //cout<<output.str();
     if (PathData.ExceededWallTime()) {
        cout << "PIMC++ exceeded wall clock limit.  Exitting LoopClass.\n";
        return;
      }
      else {
        struct timeval start, end;
        struct timezone tz;
        gettimeofday(&start, &tz);

        (*iter)->DoEvent();

        gettimeofday(&end,   &tz);
        (*iter)->TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
      }
    }
  }
}


MoveClass*
LoopClass::FindMove(string name)
{
  list<MoveClass*>::iterator iter;
  for (iter=Moves.begin(); iter != Moves.end(); iter++)
    if ((*iter)->Name == name)
      return (*iter);
  return NULL;
}

ObservableClass*
LoopClass::FindObservable(string name) 
{
  list<ObservableClass*>::iterator iter;
  for (iter=Observables.begin(); iter != Observables.end(); iter++)
    if ((*iter)->Name == name)
      return (*iter);
  return NULL;
}

void
LoopClass::Read(IOSectionClass &in)
{
  assert (in.ReadVar("Steps", NumSteps));
  if(!in.ReadVar("Equilibrate",Equilibrate))
    Equilibrate = 0;
  Read(in, NumSteps);
  ostringstream loopStream;
  loopStream << "Loop(" << NumSteps << ")";
  Name = loopStream.str();
}

void
LoopClass::Read(IOSectionClass &in, int steps)
{
  NumSteps = steps;
  int numSections = in.CountSections();
  for (int secNum=0; secNum<numSections; secNum++) {
    in.OpenSection(secNum);
    if ((in.GetName()=="Move")){
      string name;
      assert (in.ReadVar("Name", name));
      EventClass *event = FindMove(name);
      if (event == NULL) {
        cerr << "Unknown move """ << name << """.\n";
        abort();
      }
      Events.push_back(event);
    }
    else if(in.GetName() == "Switch") {
      Events.push_back(new SwitchClass(PathData, IOSection));
    }
    else if(in.GetName() == "Observe") {
      string name;
      assert (in.ReadVar("Name", name));
      EventClass *event = FindObservable(name);
      if (event == NULL) {
        cerr << "Unknown observable """ << name << """.\n";
        abort();
      }
      Events.push_back(event);
    }
    else if (in.GetName() == "Loop") {
      LoopClass *newLoop = new LoopClass(PathData, IOSection, Moves, Observables);
      newLoop->Read(in);
      Events.push_back(newLoop);
    }
    else if (in.GetName() == "WriteData")
      Events.push_back(new WriteDataClass(PathData,IOSection,Moves,Observables));
    in.CloseSection();
  }
}
