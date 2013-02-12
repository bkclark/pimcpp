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

#include "InputFile.h"

int IsLetter(char ch)
{
  return (((ch >= 'a') && (ch <= 'z')) ||
	  ((ch >= 'A') && (ch <= 'Z')));
}

int IsNumber (char ch)
{
  return ((ch >= '0') && (ch <= '9'));
}

int IsAllowed (char ch)
{
  return (IsLetter(ch) || IsNumber(ch) || (ch == '_'));
}


int IsWhiteSpace (char ch)
{
  return ((ch == ' ') || (ch == '\t') || (ch == '\n'));
}


void
InputBuffer::Write (FILE *fout)
{
  for (int i=0; i<size; i++)
    fputc ((int) buffer(i), fout);
}



//////////////////////////////////////////////////////////////
// Read a file into the buffer, ignoring C++ style comments.//
//////////////////////////////////////////////////////////////
int
InputBuffer::Read(char *FileName)
{
  int InQuotes=0;

  FILE *fin;
  if ((fin = fopen (FileName, "r")) == NULL)
    return (0);

  int count = 0;

  char ch;
  do
    {
      ch = fgetc(fin);
      if (ch != EOF)
	{
	  if (InQuotes)
	    {
	      if (ch == '\"')
		InQuotes = 0;
	      
	      count++;
	    }
	  else if (ch == '\"')
	    {
	      InQuotes = 1;
	      count++;
	    }
	  else if (ch == '/')
	    {
	      char ch2 = fgetc(fin);
	      if (ch2 == '/')
		{
		  char ch3;
		  do
		    {
		      ch3 = fgetc(fin);
		    } while ((ch3 != EOF) && (ch3 != '\n'));
		}

	      else if (ch2 == '*')
		{
		  int found = 0;
		  char ch3;
		  ch3 = fgetc(fin);
		  while (!found)
		    {
		      if (ch3 == '*')
			{
			  ch3 = fgetc(fin);
			  if (ch3 == '/')
			    found = 1;
			}
		      else
			ch3 = fgetc(fin);
		    }
		}
	      else
		{
		  count ++;
		  if (ch2 != EOF)
		    count++;

		}
	    }
	  else if (ch == '#')
	    {
	      char str[500];
	      fgets(str, 500, fin);
	      int len = strlen (str);
	      if (!strncmp(str, "include", 7))
	      {
		int i = 7;
		while ((i<len) && ((str[i]==' ') || (str[i]=='\t')))
		  i++;
		if (str[i] != '\"')
		  {
		    cerr << "Badly formed #include directive " 
			 << "in InputBuffer::Read().\n";
		    exit(1);
		  }
		else
		  {
		    i++;
		    char IncludeName[500];
		    int j=0;
		    while ((str[i] != '\"') && (i < len))
		      {
			IncludeName[j] = str[i];
			i++;
			j++;
		      }
		    IncludeName[j] = '\0';
		    if (str[i] != '\"')
		      {
			cerr << "Badly formed #include directive " 
			     << "in InputBuffer::Read().\n";
			exit(1);
		      }
		    InputBuffer IncludeBuf;
		    if (!IncludeBuf.Read(IncludeName))
		      {
			cerr << "Error reading included file "
			     << IncludeName << " in InputBuffer::Read().\n";
			exit(1);
		      }
		    //IncludeBuf.Write(stderr);
		    count += IncludeBuf.Size();
		  }
	      }
	    }
	  else
	    count++;
	}
    } while (ch != EOF);
			 
  // Now resize the buffer
  Resize (count);
  
  // Now actually read it in;
  fseek(fin, (long)0, SEEK_SET);
  count = 0;
  
  do
    {
      ch = fgetc(fin);
      if (ch != EOF)
	{
	  if (InQuotes)
	    {
	      if (ch == '\"')
		InQuotes = 0;
	      
	      buffer(count) = ch;
	      count++;
	    }
	  else if (ch == '\"')
	    {
	      InQuotes = 1;
	      buffer(count) = ch;
	      count++;
	    }
	  else if (ch == '/')
	    {
	      char ch2 = fgetc(fin);
	      if (ch2 == '/')
		{
		  char ch3;
		  do
		    {
		      ch3 = fgetc(fin);
		    } while ((ch3 != EOF) && (ch3 != '\n'));
		}
	      else if (ch2 == '*')
		{
		  int found = 0;
		  char ch3;
		  ch3 = fgetc(fin);
		  while (!found)
		    {
		      if (ch3 == '*')
			{
			  ch3 = fgetc(fin);
			  if (ch3 == '/')
			    found = 1;
			}
		      else
			ch3 = fgetc(fin);
		    }
		}
	      else
		{
		  buffer(count) = ch;
		  count ++;
		  if (ch2 != EOF)
		    {
		      buffer(count) = ch2;
		      count++;
		    }
		}
	    }
	  else if (ch == '#')
	    {
	      char str[500];
	      fgets(str, 500, fin);
	      int len = strlen (str);
	      if (!strncmp(str, "include", 7))
	      {
		int i = 7;
		while ((i<len) && ((str[i]==' ') || (str[i]=='\t')))
		  i++;
		if (str[i] != '\"')
		  {
		    cerr << "Badly formed #include directive " 
			 << "in InputBuffer::Read().\n";
		    exit(1);
		  }
		else
		  {
		    i++;
		    char IncludeName[500];
		    int j=0;
		    while ((str[i] != '\"') && (i < len))
		      {
			IncludeName[j] = str[i];
			i++;
			j++;
		      }
		    IncludeName[j] = '\0';
		    if (str[i] != '\"')
		      {
			cerr << "Badly formed #include directive " 
			     << "in InputBuffer::Read().\n";
			exit(1);
		      }
		    InputBuffer IncludeBuf;
		    IncludeBuf.Read(IncludeName);
		    for (int k=0; k<IncludeBuf.Size(); k++)
		      {
			buffer(count) = IncludeBuf(k);
			count++;
		      }
		  }
	      }

	    }
	  else
	    {
	      buffer(count) = ch;
	      count++;
	    }
	}
    } while (ch != EOF);
  fclose (fin);
  return (1);
}		      
	
//////////////////////////////////////////////////////////////
// Finds a block delimited by a starting character and      //
// ending character, eg. '{' and '}'.  If the find is       //
// successful, BlockBuffer is a buffer containing the block //
// delimited by those characters.  On success, the position //
// in the original buffer will be just after the block and  //
// a value of 1 is returned.  If it is not successful,      //
// a -1 is returned if there is a start but no end, and 0   //
// is returned if there was no start.  Note that one can    //
// have nested blocks.  eg: {hello {world}}.                //
//////////////////////////////////////////////////////////////      
int
InputBuffer::FindBlock (char StartChar, char EndChar,
			InputBuffer &BlockBuffer)
{
  int SavePos = pos;

  // First find StartChar
  while ((pos < size) && (buffer(pos) != StartChar) && 
	 (buffer(pos) != EndChar))
    pos++;
  if ((pos == size) || (buffer(pos) != StartChar))
    {
      pos = SavePos;
      return (0);
    }

  pos++;
  int StartPos = pos;

  // Now check for nested blocks
  InputBuffer TempBlock;
  int retval;
  while ((retval = FindBlock (StartChar, EndChar, TempBlock)) == 1);
  if (retval == -1)
    {
      pos = SavePos;
      return (-1);
    }
  
  // Now find EndChar
  while ((pos < size) && (buffer(pos) != EndChar))
    pos++;
  if (buffer(pos) != EndChar)
    {
      pos = SavePos;
      return (-1);
    }
  
  int EndPos = pos-1;
  pos++;
  
  BlockBuffer.size = EndPos - StartPos + 1;
  BlockBuffer.buffer.reference(buffer(blitz::Range(StartPos,EndPos)));
  BlockBuffer.pos = 0;  
  return (1);
}



int 
InputBuffer::FindQuoteBlock (InputBuffer &BlockBuffer)
{
  int SavePos = pos;
  while ((pos < size) && (buffer(pos) != '\"'))
    pos++;
  if ((pos == size) || (buffer(pos) != '\"'))
    {
      pos = SavePos;
      return (0);
    }
  
  pos++;
  int StartPos = pos;

  while ((pos < size) && (buffer(pos) != '\"'))
    pos++;
  if (buffer(pos) != '\"')
    {
      pos = SavePos;
      return (0);
    }
  int EndPos = pos-1;
  pos++;
  
  BlockBuffer.size = EndPos - StartPos + 1;
  BlockBuffer.buffer.reference(buffer(blitz::Range(StartPos,EndPos)));
  BlockBuffer.pos = 0;  
  return (1);

}


//////////////////////////////////////////////////////////////
// Return 1 is the buffer starts with str.  Returns 0       //
// otherwise.
//////////////////////////////////////////////////////////////
int
InputBuffer::StartString (char *str)
{
  int len = strlen (str);
  int match = 1;
  int i = 0;
  while ((str[i]==buffer(pos+i)) && (i<len) && ((i+pos) < size))
    i++;
  if (i<len)
    return (0); 
  return (1);
}


int
InputBuffer::FindName (char *Name)
{
   int SavePos = pos;
  
  int found = 0;
  int done = 0;
  int len = strlen(Name);
  int InQuotes = 0;
  while (!done && (pos<size))
    {
      // Look for Name
      if (buffer(pos) == '\"')
	{
	  InQuotes = !InQuotes;
	  pos++;
	}
      else if (InQuotes)
	pos++;
      else if (!InQuotes)
	if (StartString (Name))
	  {
	    if ((pos == 0) || (!IsAllowed(buffer(pos-1))))
	      if (((pos+len)<size)&&
		   (!IsAllowed(buffer(pos+len))))
		{
		  found=1;
		  done=1;
		}
	      else
		pos++;
	    else
	      pos++;
	  }
	else
	  pos++;
      if (pos == size)
	done = 1;
    }
  
  if (found)
    {
      pos += len;
      return (1);
    }
  else
    {
      pos = SavePos;
      return (0);
    }
}

int InputBuffer::FindSection (char *SectionName, InputBuffer &SectionBuffer,
			      bool rewind)
{
  int SavePos = pos;
  int found = FindName(SectionName);
  if (found)
    found = FindBlock('{', '}', SectionBuffer);
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return (found);
}
	      

int InputBuffer::FindVarBuf (char *VarName, InputBuffer &VarBuffer)
{
  int SavePos = pos;
  int found = FindName(VarName);

  if (!found)
    pos = SavePos;
  else
    {
      found = FindBlock('=', ';', VarBuffer);
      if (!found)
	pos = SavePos;
    }
  return(found);
}	       		  





int InputBuffer::ReadDouble(double &num)
{
  char tempbuf[size-pos+1];
  char *endptr;

  // Copy buffer from pos into a properly NULL terminated string;
  for (int i=pos; i<size; i++)
    tempbuf[i-pos] = buffer(i);
  tempbuf[size-pos] = '\0';

  num = strtod (tempbuf, &endptr);
  if (tempbuf == endptr)
    return (0);
  else
    {
      pos += endptr-tempbuf;
      return(1);
    }
}


int InputBuffer::ReadInt(int &num)
{
  char tempbuf[size-pos+1];
  char *endptr;

  // Copy buffer from pos into a properly NULL terminated string;
  for (int i=pos; i<size; i++)
    tempbuf[i-pos] = buffer(i);
  tempbuf[size-pos] = '\0';

  num = strtol (tempbuf, &endptr, 10);
  if (tempbuf == endptr)
    return (0);
  else
    {
      pos += endptr-tempbuf;
      return(1);
    }
}


bool IsWhite (char c)
{
  return ((c == ' ') || (c == '\t') || (c == '\n'));
}


int InputBuffer::ReadBool(bool &IsTrue)
{
  char tempbuf[size-pos+1];
  char *endptr;

  // Copy buffer from pos into a properly NULL terminated string;
  for (int i=pos; i<size; i++)
    tempbuf[i-pos] = buffer(i);
  tempbuf[size-pos] = '\0';

  // Remove initial whitespace characters;
  while (IsWhite (tempbuf[0]))
    {
      int i = 0;
      while (tempbuf[i] != '\0')
	{
	  tempbuf[i] = tempbuf[i+1];
	  i++;
	}
    }

  // Remove trailing whitespace characters
  int i = strlen (tempbuf) - 1;
  while (IsWhite (tempbuf[i]) && (i>0))    
    tempbuf[i--] = '\0';

  // make lowercase
  for (int i=0; i<strlen(tempbuf); i++)
    tempbuf[i] = tolower(tempbuf[i]);

  if (!strcmp(tempbuf, "false"))
    IsTrue = false;
  else if (!strcmp(tempbuf, "true"))
    IsTrue = true;
  else
    return 0;
  return (1);
}



/*int InputBuffer::ReadVector(blitz::Array<scalar,1> &vec)
{
  int count = 0;
  int SavePos = pos;
 
  double temp;
  while (ReadDouble(temp))
    count++;
  pos = SavePos;
  vec.resize(count);
  for (int i=0; i<count; i++)
    ReadDouble(vec(i));
  return (count);
}*/

// Faster version
int InputBuffer::ReadVector(blitz::Array<scalar,1> &vec)
{
  char tempbuf[size-pos+1];
  char *endptr;

  // Copy buffer from pos into a properly NULL terminated string;
  for (int i=pos; i<size; i++)
    tempbuf[i-pos] = buffer(i);
  tempbuf[size-pos] = '\0';

  int count = 0;

  char *bufhold = tempbuf;

  int done = 0;
  while (!done)
    {
      strtod (bufhold, &endptr);
      if (bufhold == endptr)
	done =1;
      else
	{
	  bufhold = endptr;
	  count++;
	}
    }

  vec.resize(count);
  bufhold = tempbuf;
  for (int i=0; i<count; i++)
    {
      vec(i) = strtod(bufhold, &endptr);
      bufhold = endptr;
    }
  
  pos += endptr - tempbuf;
  return (count);
}



/*int InputBuffer::ReadVector(blitz::Array<int,1> &vec)
{
  int count = 0;
  int SavePos = pos;
 
  int temp;
  while (ReadInt(temp))
    count++;
  pos = SavePos;
  vec.resize(count);
  for (int i=0; i<count; i++)
    ReadInt(vec(i));
  return (count);
}*/


int InputBuffer::ReadQuoteString (string &str)
{
  InputBuffer tempbuf;
  int found;

  found = FindQuoteBlock (tempbuf);
  str = "";
  if (found)
    for (int i=0; i<tempbuf.size; i++)
      str += tempbuf(i);
  return (found);
}


int InputBuffer::ReadVector(blitz::Array<string,1> &vec)
{
  int count = 0;
  int SavePos = pos;
 
  InputBuffer tempbuf;
  while (FindQuoteBlock(tempbuf))
    count++;
  pos = SavePos;
  vec.resize(count);
  for (int i=0; i<count; i++)
    ReadQuoteString(vec(i));
  return (count);
}



// Faster version
int InputBuffer::ReadVector(blitz::Array<int,1> &vec)
{
  char tempbuf[size-pos+1];
  char *endptr;

  // Copy buffer from pos into a properly NULL terminated string;
  for (int i=pos; i<size; i++)
    tempbuf[i-pos] = buffer(i);
  tempbuf[size-pos] = '\0';

  int count = 0;

  char *bufhold = tempbuf;

  int done = 0;
  while (!done)
    {
      strtol (bufhold, &endptr, 10);
      if (bufhold == endptr)
	done =1;
      else
	{
	  bufhold = endptr;
	  count++;
	}
    }

  vec.resize(count);
  bufhold = tempbuf;
  for (int i=0; i<count; i++)
    {
      vec(i) = strtol(bufhold, &endptr, 10);
      bufhold = endptr;
    }
  
  pos += endptr - tempbuf;
  return (count);
}






int InputBuffer::ReadString(char str[], int max)
{
  int count=0;
  while ((pos<size) && (count<(max-1)))
  {
    str[count] = buffer(pos);
    count++;
    pos++;
  }
  if(count<max)
    str[count] = '\0';
  return (count);
}


//////////////////////////////////////////////////////////////
// Named variable reading routines:                         //
//////////////////////////////////////////////////////////////

////////////
// double //
////////////
int InputBuffer::ReadVar (char *Name, double &Val, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);

  if (found)
    found = VarBuf.ReadDouble (Val);
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);
}


/////////
// int //
/////////
int InputBuffer::ReadVar (char *Name, int &Val, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);
  if (!found)
    pos = SavePos;
  else
    {
      found = VarBuf.ReadInt (Val);
      if (!found)
	pos = SavePos;
    }
  if (rewind)
    Rewind();
  return(found);
}



//////////
// bool //
//////////
int InputBuffer::ReadVar (char *Name, bool &IsTrue, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);
  if (found)
    {
      found = VarBuf.ReadBool (IsTrue);
      if (!found)
	pos = SavePos;
    }
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);
}


//////////////////////
// blitz::Array of doubles //
//////////////////////
int InputBuffer::ReadVar (char *Name, blitz::Array<double,1> &Vec, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);
  if (found)
    {
      InputBuffer VecBuf;
      found = VarBuf.FindBlock ('[', ']', VecBuf);
      if (found)
	found = VecBuf.ReadVector (Vec);
    }
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);
}

///////////////////
// blitz::Array of ints //
///////////////////
int InputBuffer::ReadVar (char *Name, blitz::Array<int,1> &Vec, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);

  if (found)
    {
      InputBuffer VecBuf;
      found = VarBuf.FindBlock ('[', ']', VecBuf);
      if (found)
	found = VecBuf.ReadVector (Vec);
    }
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);
}


//////////////////////
// Character string //
//////////////////////
int InputBuffer::ReadVar (char *Name, char *str, int maxlength, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf;
  int found = FindVarBuf (Name, VarBuf);
  if (found)
    {
      InputBuffer StrBuf;
      found = VarBuf.FindQuoteBlock (StrBuf);
      if (found)
	found = StrBuf.ReadString (str, maxlength);
    }
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);

}


////////////////////////////////
// blitz::Array of character strings //
////////////////////////////////
int InputBuffer::ReadVar (char *Name, blitz::Array<string,1> &Vec, bool rewind)
{
  int SavePos = pos;
  InputBuffer VarBuf, StrBuf;
  int found = FindVarBuf (Name, VarBuf);
  if (found)
    {
      InputBuffer VecBuf;
      found = VarBuf.FindBlock ('[', ']', VecBuf);
      if (found)
	found = VecBuf.ReadVector (Vec);
    }
  if (!found)
    pos = SavePos;
  if (rewind)
    Rewind();
  return(found);
}



//////////////////////////////////////////////////////////////
//                           Tests                          //
//////////////////////////////////////////////////////////////



void TestRead()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  InBuf.Write(stderr);
}


void TestBlock()
{
  InputBuffer InBuf;
  InputBuffer BlockBuf;
  InputBuffer Block2;
  InBuf.Read ("test.txt");
  int retval = InBuf.FindBlock('[', ']', BlockBuf);
  int retval2 = InBuf.FindBlock('[', ']', Block2);
  cerr << "Found = " << retval2 << "\n";
  Block2.Write(stderr);
}

void TestSection()
{
  InputBuffer InBuf;
  InputBuffer SectionBuf, SectionBuf2, SB3;
  InBuf.Read ("test.txt");
  
  int retval = InBuf.FindSection("level1", SectionBuf);
  cerr << "Section1:\n"; 
  SectionBuf.Write(stderr);
  retval = SectionBuf.FindSection("level2", SectionBuf2);
  cerr << "Section1:\n"; 
  SectionBuf.Write(stderr);
  
  retval = InBuf.FindSection("level1", SB3);
  cerr << "Found = " << retval << "\n";

  cerr << "Section2:\n";
  SectionBuf2.Write(stderr);

  cerr << "SB31:\n";
  SB3.Write(stderr);
}
  

void TestDouble()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  InputBuffer NumSec;

  InBuf.FindSection("Nums", NumSec);

  cerr << "NumSec = \n\"";
  NumSec.Write(stderr);
  cerr << "\"\n";


  double f, g;

  int ffound = NumSec.ReadVar("f", f);
  NumSec.Rewind();
  int gfound = NumSec.ReadVar("g", g);

  cerr << "ffound = " << ffound << "\n";
  cerr << "gfound = " << gfound << "\n";
  cerr << "f = " << f << "\n";
  cerr << "g = " << g << "\n";
}

void TestInt()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  InputBuffer NumSec;

  InBuf.FindSection("Nums", NumSec);

  cerr << "NumSec = \n\"";
  NumSec.Write(stderr);
  cerr << "\"\n";


  int f, g;

  int ffound = NumSec.ReadVar("f", f);
  NumSec.Rewind();
  int gfound = NumSec.ReadVar("g", g);

  cerr << "ffound = " << ffound << "\n";
  cerr << "gfound = " << gfound << "\n";
  cerr << "f = " << f << "\n";
  cerr << "g = " << g << "\n";
} 



void TestBool()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  InputBuffer BoolSec;

  InBuf.FindSection("Bools", BoolSec);

  cerr << "BoolSec = \n\"";
  BoolSec.Write(stderr);
  cerr << "\"\n";


  bool IsGoat, IsPig;

  int GoatFound = BoolSec.ReadVar("IsGoat", IsGoat);
  BoolSec.Rewind();
  int PigFound = BoolSec.ReadVar("IsPig", IsPig);

  cerr << "GoatFound = " << GoatFound << "\n";
  cerr << "PigFound = " << PigFound << "\n";
  cerr << "IsGoat = " << IsGoat << "\n";
  cerr << "IsPig = " << IsPig << "\n";
} 



void TestVector()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  InputBuffer NumSec;

  InBuf.FindSection("Nums", NumSec);
  cerr << "NumSec = \n\"";
  NumSec.Write(stderr);
  cerr << "\"\n";

  blitz::Array<int, 1> a;
  blitz::Array<scalar, 1> v;

  NumSec.ReadVar ("a", a);
  int found = NumSec.ReadVar("v", v);
  NumSec.Rewind();

  cerr << "found = " << found << "\n";
  cerr << "a = " << a << "\n";
  cerr << "v = " << v << "\n";
}   



void TestInclude()
{
   InputBuffer InBuf;
   InBuf.Read ("test.txt");

   InBuf.Write(stderr);
}


void TestString()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  InputBuffer NumSec;

  InBuf.FindSection("Nums", NumSec);

  char b[50];

  int found = NumSec.ReadVar("b", b, 50);
  NumSec.Rewind();

  cerr << "found = " << found << "\n";
  cerr << "b = " << b << "\n";
}   



void TestStringVector()
{
  InputBuffer InBuf;
  InBuf.Read ("test.txt");
  
  blitz::Array<string,1> stringvec;
  InputBuffer Sec;

  InBuf.FindSection("StringVector", Sec);

  int found = Sec.ReadVar("strs", stringvec);
  Sec.Rewind();

  cerr << "found = " << found << "\n";
  cerr << "stringvec = " << stringvec << "\n";
}   

#ifdef INPUT_FILE_TEST
main()
{
  TestStringVector();
  TestBool();
  TestString();
  TestVector();
  TestInt();
  TestDouble();
   TestSection();
  TestBlock();
  TestRead();
  //  TestToken();
}
#endif



