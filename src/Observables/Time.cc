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

#include "Time.h"
#include <sys/time.h>


void MCTimeClass::Accumulate()
{

  
}


void MCTimeClass::WriteBlock()
{
  if (FirstTime){
    FirstTime = false;
    BlockNumber = 0;
    MoveTime.resize(Moves.size());
    ObservableTime.resize(Observables.size());
    ActionTime.resize(Actions.size());
    MoveTime = 0;
    ObservableTime = 0;
    ActionTime = 0;
    TotalTime = 0;
    Array<string,1> moveNames (Moves.size());
    Array<string,1> observableNames (Observables.size());
    Array<string,1> actionNames (Actions.size());
    int i;
    list<MoveClass*>::iterator moveIter;
    for (moveIter=Moves.begin(),i=0; moveIter!=Moves.end(); moveIter++,i++)
      moveNames(i) = ((*moveIter)->Name);
    list<ObservableClass*>::iterator observableIter;
    for (observableIter=Observables.begin(),i=0; observableIter!=Observables.end(); observableIter++,i++)
      observableNames(i) = ((*observableIter)->Name);
    list<ActionBaseClass*>::iterator actionIter;
    for (actionIter=Actions.begin(),i=0; actionIter!=Actions.end(); actionIter++,i++)
      actionNames(i) = ((*actionIter)->GetName());

    if (PathData.Path.Communicator.MyProc()==0) {
      IOSection.WriteVar("MoveNames", moveNames);
      IOSection.WriteVar("ObservableNames", observableNames);
      IOSection.WriteVar("ActionNames", actionNames);
    }
  }
  gettimeofday(&end, &tz);
  double BlockTime = ((double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec));
  TotalTime += BlockTime;
  gettimeofday(&start, &tz);
  //  TotalTime+=(double)(clock()-StartTime)/(double)CLOCKS_PER_SEC;
  //  StartTime=clock(); << TotalTime << endl;
  TotalTimeVar.Write(TotalTime);
  list<MoveClass*>::iterator moveIter;
  int i=0;
  for (moveIter=Moves.begin();moveIter!=Moves.end();moveIter++){
    MoveTime(i)=((*moveIter)->TimeSpent)/TotalTime;
    i++;
  }
  i = 0;
  list<ObservableClass*>::iterator observeIter; 
  for (observeIter=Observables.begin();
       observeIter!=Observables.end();observeIter++){
    ObservableTime(i)=((*observeIter)->TimeSpent)/TotalTime;
    i++;
  }

  list<ActionBaseClass*>::iterator actionIter; 
  for (actionIter=Actions.begin(),i=0;
       actionIter!=Actions.end();actionIter++,i++){
    ActionTime(i)=((*actionIter)->TimeSpent)/TotalTime;
  }
  


  MoveTimeVar.Write(MoveTime);
  ObservableTimeVar.Write(ObservableTime);
  ActionTimeVar.Write(ActionTime);

  if (PathData.Path.Communicator.MyProc()==0) {
    IOSection.FlushFile();
    // Block Info
    cout << PathData.Path.CloneStr << " #: " << BlockNumber << " t: " << BlockTime << " T: " << TotalTime << endl;
  }
  BlockNumber++;
}

void MCTimeClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
  gettimeofday(&start, &tz);

}



void MCTimeClass::WriteInfo()
{
//   Array<string,1> movesString;
//   Array<string,1> observeString;
//   movesString.resize(Moves.size());
//   observeString.resize(Observables.size());
//   for (int i=0;i<Moves.size();i++){
//     movesString(i)=Moves(i)->Name;
//   }
//   for (int i=0;i<Observables.size();i++){
//     observeString(i)=Observables(i)->Name;
//   }
  //  IOSection.WriteVar("Move Names",movesString);
  //  IOSection.WriteVar("Observable Names",observeString);
  
}
