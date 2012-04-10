/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
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
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "PathDataClass.h"
#include "WriteData.h"
#include "Moves/MoveBase.h"
#include "Observables/ObservableBase.h"

void
WriteDataClass::DoEvent()
{
  list<MoveClass*>::iterator moveIter;
  for (moveIter = Moves.begin(); moveIter != Moves.end(); moveIter++) 
    (*moveIter) -> WriteRatio();

  list<ObservableClass*>::iterator observeIter;
  for (observeIter = Observables.begin(); observeIter != Observables.end(); observeIter++) 
    (*observeIter) -> WriteBlock();
  if (PathData.Path.Communicator.MyProc() ==  0)
    IOSection.FlushFile();
}


void
WriteDataClass::Read(IOSectionClass &in)
{
  //do nothing for now
}
