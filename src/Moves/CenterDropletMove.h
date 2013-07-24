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

#ifndef CENTERDROPLET_MOVE_H
#define CENTERDROPLET_MOVE_H

#include "MoveBase.h"


class CenterDropletClass : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  int Species;

  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="CenterDroplet");
      assert(moveInput.ReadVar("name",Name));
      string speciesString;
      assert(moveInput.ReadVar("species",speciesString));
      Species=PathData.Path.SpeciesNum(speciesString);

    }
  CenterDropletClass(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(PathData.Path.NumParticles());
      for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
	ActiveParticles(ptcl)=ptcl;
      /* do nothing for now */
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};



#endif
