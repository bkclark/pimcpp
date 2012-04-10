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

#include "MirroredClass.h"

ModeType ActiveCopy=OLDMODE;

void MirroredClassTest()
{
  MirroredClass<double> a;
  double b;

  b = 1.0;
  a = b;
  cerr << "a = " << a << endl;;

  Mirrored1DClass<double> c;
  Array<double,1> d;
  c.resize(3);
  d.resize(3);
  c(0) = 1.0; 
  c(1) = 2.0; c(2) = 3.0;
  d = c;
  cerr << "c = " << c.data() << endl;


  Mirrored2DClass<double> e;
  Array<double,2> f;

  e.resize(2,3);
  e.data() = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  f = e;

}
