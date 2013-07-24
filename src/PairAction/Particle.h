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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "../IO/IO.h"

using namespace IO;

class ParticleClass
{
public:
  string Name;
  double lambda;
  double Charge;
  int Ndim;
  inline void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Name", Name);
    outSection.WriteVar ("lambda", lambda);
    outSection.WriteVar ("Charge", Charge);
    outSection.WriteVar ("Ndim", Ndim);
  }
  inline bool Read  (IOSectionClass &inSection)
  {
    bool success;
    success =  inSection.ReadVar ("Name", Name);
    success &= inSection.ReadVar ("lambda", lambda);
    success &= inSection.ReadVar ("Charge", Charge);
    success &= inSection.ReadVar ("Ndim", Ndim);
    return (success);
  }
};

#endif
