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

#ifndef READ_PATH_H
#define READ_PATH_H

#include "../PathDataClass.h"
#include "MoveBase.h"

class ReadPathClass : public MoveClass
{
public:  
  int NumConfigs;
  int currConfig;
  void MakeMove();
  void Read(IOSectionClass &input);
  IOVarBase *IOVar;
  Array<int,1> ActiveParticles;
  ReadPathClass(PathDataClass &myPathData, IOSectionClass &io) :
    MoveClass(myPathData,io)
  {
    //do nothing for now
  }
};


#endif
