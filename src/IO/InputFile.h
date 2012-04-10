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

#ifndef INPUTFILE_H
#define INPUTFILE_H
//#include "../Blitz.h"
#include <blitz/array.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
typedef double scalar;

inline void Abort (char *str)
{
  cerr << str << "\nAborting.\n";
  exit(1);
}

class InputBuffer
{
private:
  int pos;
  blitz::Array<char,1> buffer;
  int size;
  int StartString (char *str);

public:

  char operator()(int i) const
  {
    return (buffer(i));
  }

  void Resize(int newsize)
  {
    size = newsize;
    buffer.resize(newsize);
  }
  
  void Rewind()
  {
    pos = 0;
  }

  int Read (char *FileName);
  void Write (FILE *fout);
  inline int Size()
  {
    return (size);
  }
  int FindBlock (char StartChar, char EndChar, InputBuffer &BlockBuffer);
  int FindQuoteBlock (InputBuffer &BlockBuffer);
  int FindName (char *Name);
  int FindSection (char *SecName, InputBuffer &SectionBuff, bool rewind=true);
  int FindVarBuf (char *VarName, InputBuffer &SectionBuff);

  int ReadDouble (double &num);
  int ReadInt (int &num);
  int ReadBool (bool &IsTrue);
  int ReadQuoteString(string &str);
  int ReadVector(blitz::Array<double,1> &vec);
  int ReadVector(blitz::Array<int,1> &vec);
  int ReadVector(blitz::Array<string,1> &vec);
  int ReadString(char str[], int max);
  int ReadVar (char *VarName, double &num, bool rewind=true);
  int ReadVar (char *VarName, int &num, bool rewind=true);
  int ReadVar (char *VarName, blitz::Array<scalar,1> &vec, bool rewind=true);
  int ReadVar (char *VarName, blitz::Array<int,1> &vec, bool rewind=true);
  int ReadVar (char *VArName, blitz::Array<string,1> &vec, bool rewind=true);
  int ReadVar (char *VarName, char *str, int maxlength, bool rewind=true);
  int ReadVar (char *VarName, bool &IsTrue, bool rewind=true);
  
  InputBuffer()
  {
    size = 0;
    pos = 0;
  }
};



#endif

