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

#include "InputOutput.h"

#include <fstream>

/// Simply prints 3*num spaces
inline void ASCIIPrintIndent(int num)
{
  for (int counter=0;counter<num*2;counter++)
    cout<<' ';
}

/// Simply prints 3*num spaces
inline void ASCIIPrintIndent(int num,ofstream &outFile)
{
  for (int counter=0;counter<num*2;counter++)
    outFile<<' ';
}


/// Prints an indented hierarchy of sections and variable names to
/// cout. 
void IOTreeASCIIClass::PrintTree(int indentNum)
{
  ASCIIPrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    ASCIIPrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

/// Calls PrintTree(0)
void IOTreeASCIIClass::PrintTree()
{
  PrintTree(0);
}




/// Returns true if theChar is a special character that should be
/// parsed into its own token.
bool isSpecial(char theChar)
{
  return ( (theChar=='(') ||
	   (theChar==')') ||
	   (theChar=='{') ||
	   (theChar=='}') ||
	   (theChar=='[') ||
	   (theChar==']') ||
	   (theChar=='<') ||
	   (theChar=='>') ||	
	   (theChar=='=') ||
	   (theChar==';') ||
	   (theChar==','));
}
	   
/// Returns true if theChar is a space, newline, tab, or carriage return.      
bool isWhiteSpace(char theChar)
{
  return ( (theChar=='\n') ||
	   (theChar==' ' ) ||
	   (theChar=='\t') ||
	   (theChar=='\r'));
}
      

/// Returns true if theChar is a letter or underscore		      
bool isAlpha(char theChar)
{
  return ((theChar>='a' && theChar<='z') || (theChar>='A' && theChar<='Z')
	  ||theChar=='_');
}

/// Returns true if theChar is a digit
bool isDigit(char theChar)
{
  return (theChar>='0' && theChar<='9');
}

/// Returns true if theChar is the a valid character for starting a
/// number.  Includes a digit, a '.' or a '-'
bool isNumStart(char theChar)
{
  return ((isDigit(theChar)) || (theChar=='.') || (theChar=='-'));
}

/// Returns true if ch is a valid character comprising a number.
bool isNumChar (char ch)
{
  return (isDigit(ch) || (ch =='.') || (ch=='e') 
	  || (ch=='-') || (ch=='+'));
}


/// Tokenize takes an array of characters and constructs a list of
/// TokenClass objects.  Each token has a string and a line number.
/// Valid tokens are special characters: "(){}[]<>,", quoted strings,
/// words, or numbers.  White space is not significant, except in
/// separating tokens.
void Tokenize(blitz::Array<char,1> buffer, list<TokenClass>& tokenList)
{
  int pos=0;
  int lineNum=1;
  while (pos<buffer.size()){
    if (isSpecial(buffer(pos))){
      TokenClass tempToken;
      tempToken.Str+=buffer(pos);
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
      pos++;
    }
    else if (buffer(pos)=='\"'){
      TokenClass tempToken;
      tempToken.Str="\"";
      pos++;
      while (buffer(pos)!='\"'){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      pos++;
      tempToken.Str+='\"';
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (isAlpha(buffer(pos))){
      TokenClass tempToken;
      while ((isAlpha(buffer(pos)) || isDigit(buffer(pos)))){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (isWhiteSpace(buffer(pos))){
      if (buffer(pos)=='\n')
	lineNum++;
      pos++;
    }
    else if (isNumStart(buffer(pos))){
      TokenClass tempToken;
      while (isNumChar(buffer(pos))){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (buffer(pos)=='\0')
      break;
    else {
      cerr<<"There was a token we do not recognize in line "<<lineNum<<endl;
      cerr <<"The rest of the file is as follows:\n";
      while (pos<buffer.size()) {
	cerr << (int)buffer(pos);
	pos++;
      }
      exit(1);
    }
  }
}	     
	  



/// Just a shortcut to look at two characters at a time.
bool checkPair(blitz::Array<char,1> &buffer,int counter,char* toSee)
{
  if (counter+1>=buffer.size()){
    return false;
  }
  if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
    return true;
  }
  else return false;

}

/// Reads a file into a character array, removing C and C++ style
/// comments. 
bool IOTreeASCIIClass::ReadWithoutComments(string fileName,
					   blitz::Array<char,1> 
					   &buffer)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile.is_open()) 
    return false;
  blitz::Array<char,1> tmpBuffer;
  int counter=0;
  bool inQuote=false;
  char dummyChar;
  while (!infile.eof()){
    infile.get(dummyChar);    
    counter++;
  }

  tmpBuffer.resize(counter);
  buffer.resize(counter);
  counter=-1;
  infile.close();
  ifstream infile2;
  infile2.open(fileName.c_str());
  while (!infile2.eof()){
    counter++;
    infile2.get(dummyChar);
    tmpBuffer(counter)=dummyChar;    
  }
  //  cout<<tmpBuffer;
  int bufferLoc=0;
  for (int counter=0;counter<tmpBuffer.size();counter++){
    if (inQuote){
      while ( counter<tmpBuffer.size() && tmpBuffer(counter) != '\"'){
	buffer(bufferLoc)=tmpBuffer(counter);
	counter++;
	bufferLoc++;
      }
      buffer(bufferLoc)=tmpBuffer(counter);
      bufferLoc++;
      inQuote=false;
    }
    else {
      if (checkPair(tmpBuffer,counter,"//")){
	while (tmpBuffer(counter)!='\n' && counter<tmpBuffer.size()){
	  counter++;
	}
	buffer(bufferLoc)=tmpBuffer(counter); //copy the \n over
	bufferLoc++;
      }
      else if (checkPair(tmpBuffer,counter,"/*")){
	while (!checkPair(tmpBuffer,counter,"*/") && counter<tmpBuffer.size()){
	  counter++;
	  if (tmpBuffer(counter)=='\n'){
	    buffer(bufferLoc)=tmpBuffer(counter);
	    bufferLoc++;
	  }		   
	}
	counter++; //end up in the / of comment
      }
      else if (tmpBuffer(counter)=='\"'){
	inQuote=true;
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
      else {
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
    }
  }
  buffer.resizeAndPreserve(bufferLoc);
  infile2.close();
  return (true);
}



/// If isError is true, then we print out an error message giving the
/// line number and the string passed to us in ErrorStr.
inline void ReadAbort (bool isError, int lineNumber, string ErrorStr)
{
  if (isError) {
    cerr << "Error in input file at line number " << lineNumber 
	 << ":\n";
    cerr << ErrorStr;
    exit(1);
  }
}

/// Removes all double quotes from the input string and return it.
string StripQuote(string str)
{
  string newString;
  int i=0;
  assert (str[0] == '\"');
  assert (str[str.length()-1] == '\"');
  while (i<str.length())
    {
      if (str[i] != '\"')
	newString += str[i];
      i++;
    }
  return newString;
}


/// Looks at the string passed to it and returns the corresponding
/// enumerated type.  If the type is not recognized, it returns
/// NOT_ATOMIC.  
AtomicType GetType (string typeString)
{
  if (typeString=="double")
    return DOUBLE_TYPE;
  else if (typeString=="int")
    return INT_TYPE;
  else if (typeString=="string")
    return STRING_TYPE;
  else if (typeString=="bool")
    return BOOL_TYPE;
  else return NOT_ATOMIC;
}


/// Takes a token and reads its value into a double, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,double &d)
{

  char* endPtr;
  d=strtod(token.Str.c_str(),&endPtr);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Double\n");
}

/// Takes a token and reads its value into an int, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,int &d)
{

  char* endPtr;
  d=strtol(token.Str.c_str(),&endPtr,10);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Int\n");
}

/// Takes a token and reads its value into a string, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,string &d)
{

  ReadAbort (token.Str[0] != '\"', token.LineNumber, 
	     "Expected '\"'.");
  ReadAbort (token.Str[token.Str.length()-1] != '\"', token.LineNumber, 
	     "Expected '\"'.");
  d=StripQuote(token.Str);
}

/// Takes a token and reads its value into a bool, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token, bool &b)
{
  if (token.Str=="true"){
    b=true;
  }
  else if (token.Str=="false"){
    b=false;
  }
  else ReadAbort(true,token.LineNumber,"Expected true or false\n");
}

/// This template function reads a 1-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   blitz::Array<T,1> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int counter=0;counter<valArray.extent(0)-1;counter++){
    ReadAtomicVar(*iter,valArray(counter));
    iter++;
    ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
    iter++;
  }
  //Read last value
  ReadAtomicVar(*iter,valArray(valArray.extent(0)-1));
  iter++;
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


/// This template function reads a 2-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.  The data is
/// read row-ordered, i.e. the first index changes fastests as we read
/// in the values.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   blitz::Array<T,2> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int j=0; j<valArray.extent(0); j++)
    for (int i=0;i<valArray.extent(1);i++) {
      ReadAtomicVar(*iter,valArray(j,i));
      iter++;
      // Read comma if this isn't the last value.
      if ((i!=valArray.extent(1)-1) || 
	  (j!=valArray.extent(0)-1)) {
	ReadAbort(iter->Str != ",", iter->LineNumber, 
		  "Expected , not found\n");
	iter++;
      }
    }
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


/// This template function reads a 3-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.  The data is
/// read row-ordered, i.e. the first index changes fastests as we read
/// in the values.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   blitz::Array<T,3> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int k=0; k<valArray.extent(0); k++)
    for (int j=0; j<valArray.extent(1); j++)
      for (int i=0;i<valArray.extent(2);i++) {
	ReadAtomicVar(*iter,valArray(k,j,i));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(2)-1) || 
	    (j!=valArray.extent(1)-1) ||
	    (k!=valArray.extent(0)-1)) {
	  ReadAbort(iter->Str != ",", iter->LineNumber, 
		    "Expected , not found\n");
	  iter++;
	}
      }
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


/// This template function reads a 4-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.  The data is
/// read row-ordered, i.e. the first index changes fastests as we read
/// in the values.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   blitz::Array<T,4> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int k=0; k<valArray.extent(0); k++) 
    for (int j=0; j<valArray.extent(1); j++)
      for (int i=0;i<valArray.extent(2);i++) 
	for (int h=0;h<valArray.extent(3);h++) {
	  ReadAtomicVar(*iter,valArray(k,j,i,h));
	  iter++;
	  // Read comma if this isn't the last value.
	  if ((h!=valArray.extent(3)-1) || 
	      (i!=valArray.extent(2)-1) || 
	      (j!=valArray.extent(1)-1) ||
	      (k!=valArray.extent(0)-1)) {
	    ReadAbort(iter->Str != ",", iter->LineNumber, 
		      "Expected , not found\n");
	    iter++;
	  }
	}
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


VarASCIIClass *NewASCIIVar (AtomicType newType, int ndim,
			    blitz::Array<int,1> dims)
{
  if (ndim == 0) {
    if (newType == DOUBLE_TYPE)
      return new VarASCIIdouble0Class;
    else if (newType == INT_TYPE)
      return new VarASCIIint0Class;
    else if (newType == STRING_TYPE)
      return new VarASCIIstring0Class;
    else if (newType == BOOL_TYPE)
      return new VarASCIIbool0Class;
  }
  else if (ndim == 1) {
    if (newType == DOUBLE_TYPE)
      {
	VarASCIIdouble1Class *newVar = new VarASCIIdouble1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarASCIIint1Class *newVar = new VarASCIIint1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarASCIIstring1Class *newVar = new VarASCIIstring1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarASCIIbool1Class *newVar = new VarASCIIbool1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
  }
  else if (ndim == 2) {
    if (newType == DOUBLE_TYPE)
      {
	VarASCIIdouble2Class *newVar = new VarASCIIdouble2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarASCIIint2Class *newVar = new VarASCIIint2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarASCIIstring2Class *newVar = new VarASCIIstring2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarASCIIbool2Class *newVar = new VarASCIIbool2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
  }  
  else if (ndim == 3) {
    if (newType == DOUBLE_TYPE)
      {
	VarASCIIdouble3Class *newVar = new VarASCIIdouble3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarASCIIint3Class *newVar = new VarASCIIint3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarASCIIstring3Class *newVar = new VarASCIIstring3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarASCIIbool3Class *newVar = new VarASCIIbool3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
  }
  else if (ndim == 4) {
    if (newType == DOUBLE_TYPE)
      {
	VarASCIIdouble4Class *newVar = new VarASCIIdouble4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarASCIIint4Class *newVar = new VarASCIIint4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarASCIIstring4Class *newVar = new VarASCIIstring4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarASCIIbool4Class *newVar = new VarASCIIbool4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
  }  
  
}



/// Reads an array from a list of tokens, starting at the token
/// pointed to by iter.  It places the array into the newVar object.
/// It expects to begin reading after the word "Array".
VarASCIIClass * ReadArray(list<TokenClass>::iterator &iter,
			  list<TokenClass> &tokenList)
{
  ReadAbort(iter->Str != "<", iter->LineNumber, "Expected < not found\n");
  iter++;
  AtomicType myType=GetType(iter->Str);
  ReadAbort(myType==NOT_ATOMIC,iter->LineNumber,
	    "Array does not have atomic type\n");
  iter++;
  ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
  iter++;
  int numDim;
  ReadAtomicVar(*iter,numDim);
  iter++;
  ReadAbort(iter->Str != ">", iter->LineNumber, "Expected , not found\n");
  iter++;
  
  blitz::Array<int,1> dimSize(numDim);
  
  string myName=iter->Str;

  iter++;
  ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
  iter++;
  for (int counter=0;counter<numDim-1;counter++){
    ReadAtomicVar(*iter,dimSize(counter));
    iter++;
    ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
    iter++;
  }

  //Read the last dimension
  ReadAtomicVar(*iter,dimSize(numDim-1));
  iter++;
  ReadAbort(iter->Str != ")", iter->LineNumber, "Expected ) not found\n");
  iter++;
  ReadAbort(iter->Str!="=",iter->LineNumber,"Expected = not found\n");
  iter++;
  
  VarASCIIClass *newVar = NewASCIIVar (myType, numDim, dimSize);
  
  newVar->Name = myName;
  newVar->Dim=numDim;
  newVar->Type=myType;
  if (numDim == 1) {
    if (myType == DOUBLE_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIdouble1Class *)newVar)->Value);
    else if (myType == INT_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIint1Class *)newVar)->Value);
    else if (myType == STRING_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIstring1Class *)newVar)->Value);
    else if (myType == BOOL_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIbool1Class *)newVar)->Value);
  }
  if (numDim == 2) {
    if (myType == DOUBLE_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIdouble2Class *)newVar)->Value);
    else if (myType == INT_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIint2Class *)newVar)->Value);
    else if (myType == STRING_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIstring2Class *)newVar)->Value);
    else if (myType == BOOL_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIbool2Class *)newVar)->Value);
  }
  if (numDim == 3) {
    if (myType == DOUBLE_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIdouble3Class *)newVar)->Value);
    else if (myType == INT_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIint3Class *)newVar)->Value);
    else if (myType == STRING_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIstring3Class *)newVar)->Value);
    else if (myType == BOOL_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIbool3Class *)newVar)->Value);
  }
  if (numDim == 4) {
    if (myType == DOUBLE_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIdouble4Class *)newVar)->Value);
    else if (myType == INT_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIint4Class *)newVar)->Value);
    else if (myType == STRING_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIstring4Class *)newVar)->Value);
    else if (myType == BOOL_TYPE)
      ReadArrayData(iter, tokenList, ((VarASCIIbool4Class *)newVar)->Value);
  }
  return (newVar);
}



/// This function parses a variable assigment from the list of tokens,
/// creates a new VarASCIIClass object and puts the appropriate value
/// in that object.  It recognizes any of the atomic types or arrays
/// of theose atomic types.
VarASCIIClass* ReadASCIIVar (list<TokenClass>::iterator &iter,
			     list<TokenClass> &tokenList)
{
  VarASCIIClass *newVar;
  AtomicType myType=GetType(iter->Str);
  if (myType==NOT_ATOMIC){
    ReadAbort(iter->Str!="Array",iter->LineNumber,
	      "Invalid Type: "+iter->Str+"\n");
    iter++;
    newVar = ReadArray(iter,tokenList);
  }
  else {
    iter++;
    string myName=iter->Str;

    iter++;
    ReadAbort(iter->Str!="=",iter->LineNumber,"Expected equals sign\n");
    iter++;
    TokenClass valToken=*iter;
    iter++;
    ReadAbort(iter->Str!=";",iter->LineNumber,"Expected semicolon\n");
    iter++;
    blitz::Array<int,1> dims(1);
    newVar = NewASCIIVar (myType, 0, dims);
    newVar->Type=myType;
    newVar->Dim=0;
    newVar->Name = myName;
    if (myType == DOUBLE_TYPE)
      ReadAtomicVar (valToken, ((VarASCIIdouble0Class *)newVar)->Value);
    else if (myType == INT_TYPE)
      ReadAtomicVar (valToken, ((VarASCIIint0Class *)newVar)->Value);
    else if (myType == STRING_TYPE)
      ReadAtomicVar (valToken, ((VarASCIIstring0Class *)newVar)->Value);
    else if (myType == BOOL_TYPE)
      ReadAtomicVar (valToken, ((VarASCIIbool0Class *)newVar)->Value);
  }
  return(newVar);
}


/// ReadSection parses a section in the input file.  It takes as
/// arguments this sections parent, its name, a tokenlist iterator,
/// the tokenlist, and a bool telling us if we want to look for a "}"
/// at the end of the input.  If we don't, we keep parsing until the
/// buffer runs out.  Calls itself recursively as necessary, builing
/// up a tree of sections and variables.
bool IOTreeASCIIClass::ReadSection (IOTreeClass *parent,
				       string myName,
				       list<TokenClass>::iterator &iter,
				       list<TokenClass> &tokenList,
				       bool wantEndBrace)
{
  Parent = parent;
  Name = myName;
  while ((iter != tokenList.end()) && (iter->Str != "}")) {
    if (iter->Str == "Section") {
      IOTreeClass *newTree;
      iter++;
      ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
      iter++;
      string newName = iter->Str;
      iter++;
      // Check for included section
      if (iter->Str == ",") {
	// Get filename
	iter++;
	string fileName = StripQuote(iter->Str);
	iter++;
	ReadAbort (iter->Str!=")", iter->LineNumber, "Expected ) not found\n");
	iter++;
	ReadAbort (iter->Str!=";", iter->LineNumber, "Expected ; not found\n");
	iter++;	
	newTree = ReadTree (fileName, newName, this);
      }
      else {
	ReadAbort(iter->Str != ")", iter->LineNumber, 
		  "Expected ) not found\n");
	iter++;
	ReadAbort(iter->Str != "{", iter->LineNumber, 
		  "Expected { not found\n");
	iter++;
	newTree = new IOTreeASCIIClass();
	((IOTreeASCIIClass*)newTree)->ReadSection((IOTreeClass*)this,
						     newName,iter,tokenList,
						     true);         
      }
      SectionList.push_back(newTree);
    }
    else {
      VarClass *newVar =  ReadASCIIVar(iter, tokenList);
      VarList.push_back(newVar);
    }
  }
  if ((iter==tokenList.end()) && wantEndBrace) {
    cerr << "Unexpected end of file before } \n";
    exit (1);
  }
	    
  if (iter!=tokenList.end())  
    iter++;
  return (true);
}

IOTreeClass* IOTreeASCIIClass::NewSection(string name)
{
  IOTreeClass* tempSection=new IOTreeASCIIClass();
  tempSection->Name=name;
  tempSection->Parent=this;
  tempSection->MyNumber=CurrSecNum;
  CurrSecNum++;
  SectionList.push_back(tempSection);
  MarkModified();
  return tempSection;
}

void IOTreeASCIIClass::IncludeSection(IOTreeClass *newSection)
{
  newSection->MyNumber=CurrSecNum++;
  SectionList.push_back(newSection);
  MarkModified();
}



bool IOTreeASCIIClass::NewFile (string fileName,
				string mySectionName,
				IOTreeClass *parent)
{
  FileName=fileName;
  Parent=parent;
  Name=mySectionName;
  return true;
}



/// OpenFile takes a filename to open, the name of this section and
/// the parent of this section.  It reads the file into a buffer,
/// converts it to a list of tokens, then parses the tokens,
/// constructing a tree of sections containing variables lists.  
bool IOTreeASCIIClass::OpenFile(string fileName, string myName, 
				IOTreeClass *parent)
{
  blitz::Array<char,1> buffer;
  bool success = ReadWithoutComments(fileName,buffer);
  if (!success)
    return false;
  list<TokenClass> tokenList;
  Tokenize(buffer,tokenList);
  list<TokenClass>::iterator iter=tokenList.begin();
  ReadSection(parent,myName,iter,tokenList, false);
  return true;
}



void IOTreeASCIIClass::WriteSection(ofstream &outFile,int indentNum)
{
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    ASCIIPrintIndent(indentNum,outFile);
    ((VarASCIIClass*)(*varIter))->Print(outFile);    
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    if ((*secIter)->FileName==""){
      ASCIIPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<")\n";
      ASCIIPrintIndent(indentNum,outFile);
      outFile<<"{\n";
      ((IOTreeASCIIClass*)(*secIter))->WriteSection(outFile,indentNum+1);
      ASCIIPrintIndent(indentNum,outFile);
      outFile<<"}\n\n";
    }
    else {
      ASCIIPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<", \"";
      outFile<<(*secIter)->FileName<<"\");"<<endl;
      (*secIter)->FlushFile();
    }
    secIter++;
  }
}



void IOTreeASCIIClass::FlushFile()
{
  ofstream outfile;
  if ((FileName!="") && IsModified){
    outfile.open(FileName.c_str());
    WriteSection(outfile,0);
  }

  list<IOTreeClass*>::iterator iter = SectionList.begin();
  while (iter != SectionList.end()) {
    (*iter)->FlushFile();
    iter++;
  }
}

/// CloseFile recursively destroys the tree of data we constructed.
void IOTreeASCIIClass::CloseFile()
{  
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
  
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty())
    {
      SectionList.front()->CloseFile();
      delete SectionList.front();
      SectionList.pop_front();
    }
}



//////////////////////////////////////////////////////////////////////
//                             ReadInto's                           //
//////////////////////////////////////////////////////////////////////

bool VarASCIIClass::ReadInto (double &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (int &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (string &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (bool &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<double,1> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<double,2> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<double,3> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<double,4> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<int,1> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<int,2> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<int,3> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<int,4> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<string,1> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<string,2> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<string,3> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<string,4> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<bool,1> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<bool,2> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<bool,3> &val)
{ ComplainReadInto(); return false; }
bool VarASCIIClass::ReadInto (blitz::Array<bool,4> &val)
{ ComplainReadInto(); return false; }

bool VarASCIIClass::Append (double val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (int val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (string val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (bool val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<double,1> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<double,2> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<double,3> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<int,1> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<int,2> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<int,3> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<string,1> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<string,2> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<string,3> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<bool,1> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<bool,2> &val)
{ ComplainAppend(); return false; }
bool VarASCIIClass::Append (blitz::Array<bool,3> &val)
{ ComplainAppend(); return false; }


bool VarASCIIdouble0Class::ReadInto (double &val)
{ val = Value; return true;}
bool VarASCIIint0Class::ReadInto (int &val)
{  val = Value; return true; }
bool VarASCIIstring0Class::ReadInto (string &val)
{  val = Value; return true; }
bool VarASCIIbool0Class::ReadInto (bool &val)
{ val = Value; return true; }
bool VarASCIIdouble1Class::ReadInto (blitz::Array<double,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarASCIIdouble2Class::ReadInto (blitz::Array<double,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarASCIIdouble3Class::ReadInto (blitz::Array<double,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); val = Value; return true; }
bool VarASCIIdouble4Class::ReadInto (blitz::Array<double,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
  val = Value; return true; }
bool VarASCIIint1Class::ReadInto (blitz::Array<int,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarASCIIint2Class::ReadInto (blitz::Array<int,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarASCIIint3Class::ReadInto (blitz::Array<int,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarASCIIint4Class::ReadInto (blitz::Array<int,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
 val = Value; return true; }
bool VarASCIIstring1Class::ReadInto (blitz::Array<string,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarASCIIstring2Class::ReadInto (blitz::Array<string,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarASCIIstring3Class::ReadInto (blitz::Array<string,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarASCIIstring4Class::ReadInto (blitz::Array<string,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
 val = Value; return true; }
bool VarASCIIbool1Class::ReadInto (blitz::Array<bool,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarASCIIbool2Class::ReadInto (blitz::Array<bool,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarASCIIbool3Class::ReadInto (blitz::Array<bool,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarASCIIbool4Class::ReadInto (blitz::Array<bool,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
 val = Value; return true; }


/*****************************************************************
 *                         Variable  Appends                     *
 *****************************************************************/
bool VarASCIIdouble1Class::Append (double val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarASCIIint1Class::Append (int val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarASCIIstring1Class::Append (string val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarASCIIbool1Class::Append (bool val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarASCIIdouble2Class::Append (blitz::Array<double,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIdouble3Class::Append (blitz::Array<double,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIdouble4Class::Append (blitz::Array<double,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIint2Class::Append (blitz::Array<int,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIint3Class::Append (blitz::Array<int,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIint4Class::Append (blitz::Array<int,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIstring2Class::Append (blitz::Array<string,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIstring3Class::Append (blitz::Array<string,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIstring4Class::Append (blitz::Array<string,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIbool2Class::Append (blitz::Array<bool,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIbool3Class::Append (blitz::Array<bool,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarASCIIbool4Class::Append (blitz::Array<bool,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == o);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}


void VarASCIIdouble0Class::Print (ofstream &outFile)
{ outFile << "double " << Name << " = " << Value << ";\n";}
void VarASCIIint0Class::Print (ofstream &outFile)
{ outFile << "int " << Name << " = " << Value << ";\n";} 
void VarASCIIstring0Class::Print (ofstream &outFile)
{ outFile << "string " << Name << " = " << "\"" << Value << "\";\n";} 
void VarASCIIbool0Class::Print (ofstream &outFile)
{ outFile << "bool " << Name <<" = "<< (Value ? "true;\n" : "false;\n"); }
void VarASCIIdouble1Class::Print (ofstream &outFile)
{ 
  outFile << "blitz::Array<double,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarASCIIdouble2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,2> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarASCIIdouble3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarASCIIdouble4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,4> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ","
	  << Value.extent(3) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	outFile << Value(i,j,k,l) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1,
		   Value.extent(3)-1) 
	  << "];\n";
}
void VarASCIIint1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarASCIIint2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarASCIIint3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarASCIIint4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,4> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) 
	  << "," << Value.extent(3)<< ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(2)); l++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	outFile << Value(i,j,k,l) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1,
		   Value.extent(3)-1)
	  << "];\n";
}
void VarASCIIstring1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << "\"" << Value(i) << "\"" << ", ";
  outFile << "\"" << Value(Value.extent(0)-1) << "\"];\n";
}
void VarASCIIstring2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << "\"" << Value(i,j) << "\", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarASCIIstring3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << "\"" << Value(i,j,k) << "\", ";
  outFile << "\"" 
	  << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "\"];\n";
}
void VarASCIIstring4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,4> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2)
	  << "," << Value.extent(3) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1))|| (l<(Value.extent(3)-1)))
	outFile << "\"" << Value(i,j,k,l) << "\", ";
  outFile << "\"" 
	  << Value(Value.extent(0)-1,Value.extent(1)-1,
		   Value.extent(2)-1,Value.extent(3)-1) 
	  << "\"];\n";
}
void VarASCIIbool1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << (Value(i) ? "true" : "false") << ", ";
  outFile << (Value(Value.extent(0)-1) ? "true" : "false") << "];\n";
}
void VarASCIIbool2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << (Value(i,j) ? "true" : "false") << ", ";
  outFile << 
    (Value(Value.extent(0)-1,Value.extent(1)-1) ? "true" : "false") << "];\n";
}
void VarASCIIbool3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << (Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) ?
	      "true" : "false") << "];\n";
}
void VarASCIIbool4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,4> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2)
	  << "," << Value.extent(3) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	outFile << Value(i,j,k,l) << ", ";
  outFile << (Value(Value.extent(0)-1,Value.extent(1)-1,
		    Value.extent(2)-1,Value.extent(3)-1) ?
	      "true" : "false") << "];\n";
}




//////////////////////////////////////////////////////////////////////
//                             Append's                             //
//////////////////////////////////////////////////////////////////////


void IOTreeASCIIClass::WriteVar(string name, double val)
{
  VarASCIIdouble0Class *newVar = new VarASCIIdouble0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<double,1> &val)
{
  VarASCIIdouble1Class *newVar = new VarASCIIdouble1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<double,2> &val)
{
  VarASCIIdouble2Class *newVar = new VarASCIIdouble2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<double,3> &val)
{
  VarASCIIdouble3Class *newVar = new VarASCIIdouble3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<double,4> &val)
{
  VarASCIIdouble4Class *newVar = new VarASCIIdouble4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeASCIIClass::WriteVar(string name, int val)
{
  VarASCIIint0Class *newVar = new VarASCIIint0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<int,1> &val)
{
  VarASCIIint1Class *newVar = new VarASCIIint1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<int,2> &val)
{
  VarASCIIint2Class *newVar = new VarASCIIint2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<int,3> &val)
{
  VarASCIIint3Class *newVar = new VarASCIIint3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<int,4> &val)
{
  VarASCIIint4Class *newVar = new VarASCIIint4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeASCIIClass::WriteVar(string name, string val)
{
  VarASCIIstring0Class *newVar = new VarASCIIstring0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<string,1> &val)
{
  VarASCIIstring1Class *newVar = new VarASCIIstring1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<string,2> &val)
{
  VarASCIIstring2Class *newVar = new VarASCIIstring2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<string,3> &val)
{
  VarASCIIstring3Class *newVar = new VarASCIIstring3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<string,4> &val)
{
  VarASCIIstring4Class *newVar = new VarASCIIstring4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
 
void IOTreeASCIIClass::WriteVar(string name, bool val)
{
  VarASCIIbool0Class *newVar = new VarASCIIbool0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<bool,1> &val)
{
  VarASCIIbool1Class *newVar = new VarASCIIbool1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<bool,2> &val)
{
  VarASCIIbool2Class *newVar = new VarASCIIbool2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<bool,3> &val)
{
  VarASCIIbool3Class *newVar = new VarASCIIbool3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeASCIIClass::WriteVar(string name, blitz::Array<bool,4> &val)
{
  VarASCIIbool4Class *newVar = new VarASCIIbool4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
