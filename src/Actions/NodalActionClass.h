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

#ifndef NODAL_ACTION_CLASS_H
#define NODAL_ACTION_CLASS_H

#include "ActionBase.h"
#include "../Splines/CubicSpline.h"

  class PathClass;

typedef enum { FREE_PARTICLE, GROUND_STATE, GROUND_STATE_FP } NodeType;

class NodalActionClass : public ActionBaseClass
{
public:
  virtual bool IsPositive (int slice) = 0;
  //  virtual double Det(int slice)       = 0;
  //  virtual Array<double,2> GetMatrix (int slice=0) = 0;
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  virtual void Init();
  virtual bool IsGroundState() = 0;
  virtual NodeType Type() = 0;
  virtual void Setk (Vec3 kVec);
  virtual void Update();
  NodalActionClass (PathDataClass &pathData) :
    ActionBaseClass (pathData)
  {

  }
};


#endif
