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

#ifndef PA_FIT_H
#define PA_FIT_H

#include "PAclassicalFit.h"
#include "PAzeroFit.h"
#include "DavidPAClass.h"

inline PairActionFitClass *ReadPAFit (IOSectionClass &in, 
				      double smallestBeta, int numLevels)
{
  //cerr<<"Now IN ReadPAFit"<<endl;
  assert (in.OpenSection("Fits"));
  string type;
  assert (in.ReadVar("Type", type));
  //cerr<<"The type is "<<type<<endl;
  in.CloseSection (); // "Fits"
  PairActionFitClass *fit;
  //cerr<<"CHECKING "<<(type=="DavidFit")<<" "<<(1==1)<<endl;
  if (type == "classical"){
    //cerr<<"Classical fit"<<endl; 
    fit = new PAclassicalFitClass;
  }
  else if (type == "zerofit"){
    //cerr<<"Zero fit"<<endl;
    fit=new PAzeroFitClass;
  }
  else if (type=="DavidFit"){
    //cerr<<"DavidPAClass fit"<<endl;
    fit = new DavidPAClass;
  }
  else {
    cerr << "Unrecognize pair action fit type \"" 
	 << type << "\".  Exitting.\n";
    exit(1);
  }
  //cerr<<"Pre-READ"<<endl;
  fit->Read(in, smallestBeta, numLevels);
  //cerr<<"EXITING"<<endl;
  return (fit);
}


#endif
