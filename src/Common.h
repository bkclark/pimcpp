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

#ifndef COMMON_H
#define COMMON_H

// #define SIMPLE_SPRNG

/// Number of dimensions that dVec uses.
#define NDIM 2

#include <sprng.h>
#include "Blitz.h"
#include <fstream>
typedef TinyVector<double,NDIM> dVec;

/// These are the different mode types for the MirroredArrayClass
//typedef enum {OLDMODE, NEWMODE, BOTHMODE} ModeType;

///ParticleID=(species,particle number)

typedef TinyVector<int,2> ParticleID;




///from codepedia.com
bool fileExists(const std::string& fileName);


/// These are the global variables to be used to decide what part of
/// the mirrored array we are writing to and reading from  
extern int Write1;
extern int Write2; 



/// Changes the mode the entire code is running in.
//void SetMode(ModeType);




class ImageNumClass
{
public:
  int ImageNum;
  inline ImageNumClass operator-() const
  {
    ImageNumClass minusNum;
    int mask = ~((~0)<<(2*NDIM));
    int sub = 0x55555555;  // binary: 01010101010101
    minusNum.ImageNum= (((~ImageNum)-sub)&mask);
    return minusNum;
  }
  inline operator int() 
  {
    return (ImageNum);
  }
  inline ImageNumClass(int i)
  {
    ImageNum = i;
  }
  inline ImageNumClass()
  { }
};


#endif
