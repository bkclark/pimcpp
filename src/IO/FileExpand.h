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

#ifndef FILE_EXPAND_H
#define FILE_EXPAND_H

#ifndef MAC
#include <wordexp.h>
#endif
#include <stdlib.h>

inline string ExpandFileName(string fname)
{
  string outName;
#ifdef MAC
  if (fname[0] == '~') {
    outName = getenv ("HOME"); 
    fname = fname.erase(0,1);
  }
  outName.append(fname);
#else
  wordexp_t words;
  wordexp (fname.c_str(), &words, 0);
  outName = words.we_wordv[0];
  wordfree(&words);
#endif
  return outName;
}

#endif
