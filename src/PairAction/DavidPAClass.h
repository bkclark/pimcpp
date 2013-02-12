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

#ifndef DavidPAClass_H
#define DavidPAClass_H

// #include "../MPI/Communication.h" //include not needed
#include "PAFitBase.h"
#include "../Blitz.h"
#include "../Splines/LinearSpline.h"
#include <fstream>
#include <iostream>
//#include "../../PathDataClass.h"


///This is the pair action class. It uses the following formula in
///order to calculate the pair action
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
class DavidPAClass : public PairActionFitClass 
{

 private:
  
  ///Holds the Ukj coefficients for a given q
  blitz::Array<double,1> TempukjArray;
  blitz::Array<double,1> TempdukjArray;
  int NumTerms;
  //  DistanceTableClass *DistanceTable;
  /// Skips to the next string in the file whose substring matches skipToString
  inline string SkipTo(ifstream &infile, string skipToString);
  /// Reads a Fortran 3 tensor
  inline void ReadFORTRAN3Tensor(ifstream &infile, blitz::Array<double,3> &tempUkj);
 public:
  double Vimage;
  double kCutoff;
  int dummyPairActionClass;
  blitz::Array<double,1> Potential; 
  string type1,type2;
  bool SamplingTableRead;
  inline bool Read(IOSectionClass &IOSection,double desiredTau, int numLevels);
  inline void Print();
  double Udiag(double q, int level);
  double dUdiag_fast(double q, int level);
  double UDiag_exact(double q,int level);
  
  double DesiredTau;
  int TauPos;
  int NumLevels;
  int NumTau;
  int LongRangeDim;
  double LongRangeMass1;
  double LongRangeMass2;
  bool HasLongRange;
  Array<double,1> LongRangeBox;
  /// This stores the coefficients of the expansion specified above.
  /// The array index is over the levels.  You call it with the q
  /// value and a temporary array to get all of the values in that
  /// column. 
  blitz::Array<MultiCubicSpline,1> ukj; ///<(level )
  blitz::Array<LinearSpline,1> UdiagSpline;
  blitz::Array<LinearSpline,1> dUdiagSpline;
  ///Same as ukj but stores the beta derivatives.
  blitz::Array<MultiCubicSpline,1> dukj; ///<(level )
  blitz::Array<CubicSpline,1> dUdRTimesSigmaSpline; ///<(level
  blitz::Array<CubicSpline,1> dUdRTimesSigmaSpline_movers; ///<(level
  blitz::Array<CubicSpline,1> d2UdR2TimesSigmaSpline; ///<(level
  blitz::Array<CubicSpline,1> d2UdR2TimesSigmaSpline_movers; ///<(level

  blitz::Array<double,1> SamplingTau;

  /////  MultiCubicSpline Pot;
  /// Calculate the U(s,q,z) value when given s,q,z and the level 
  void calcUsqz(double s,double q,double z,int level,
		       double &U, double &dU, double &V);
  void calcUsqzFast(double s,double q,double z,int level,
		    double &U);


  /// This is the order of the fit to use. 
  int n;
  /// This is the temperature 
  //////  double tau;
  /// Function to read David's squarer file input.
  void ReadLongRangeHDF5(IOSectionClass &in);
  void ReadDavidSquarerFile(string DMFile);
  void ReadDavidSquarerFileHDF5(string DMFile);
  void ReadSamplingTable(string fileName);
  void PrintVals(double begin, double end, double dx);
  double U (double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  double V(double r);
  double dUdRTimesSigma(double r,int level);
  double dUdRTimesSigma_movers(double r,int level);
  double d2UdR2TimesSigma(double r,int level);
  double d2UdR2TimesSigma_movers(double r,int level);


  blitz::Array<double,1> kVals;
  blitz::Array<double,1> uk_long;


  /// The q-derivative of the above
  double Udiag_p(double q, int level);
  /// The q-derivative of the above
  double Udiag_pp(double q, int level);
  /// The beta-derivative of the diagonal action
  double dUdiag    (double q, int level);
  /// The q-derivative of the above
  double dUdiag_p  (double q, int level);
  /// The q-derivative of the above
  double dUdiag_pp (double q, int level);
  void Derivs   (double q, double z, double s2, int level,
	         double &d_dq, double &d_dz);
  void DerivsFD (double q, double z, double s2, int level,
		 double &d_dq, double &d_dz);
  void Derivs   (double q, double z, double s2, int level,
	         double &d_dq, double &d_dz, double &d_ds);
  void DerivsFD (double q, double z, double s2, int level,
		 double &d_dq, double &d_dz, double &d_ds);
  bool IsLongRange(); 

  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  //  void Error (Rho &rho, double &Uerror, double &dUerror);
  //  void DoFit (Rho &rho);
  void WriteFit(IOSectionClass &outSection);
};


inline bool DavidPAClass::Read(IOSectionClass &in,double x, int y)
{
  Name="DavidPAClass";
  //cerr<<"In DavidPAClass Read"<<endl;
  //cerr<<"values are "<<x<<" "<<y<<endl;
  string fileName;
  DesiredTau=x;
  NumLevels=y;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.ReadVar("NumOffDiagonalTerms",NumTerms));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;

  assert(in.ReadVar("Daviddmfile",fileName));
  int l = fileName.size();
  // slice last two characters of the filename. i.e. the extension
  string extension(fileName, l-2, 2);
  if(extension == "dm")
    ReadDavidSquarerFile(fileName.c_str());
  else if (extension == "h5")
    ReadDavidSquarerFileHDF5(fileName.c_str());
  else
    assert(0);
  string samplingTableFile;
  SamplingTableRead=in.ReadVar("SamplingTableFile",samplingTableFile);
  if (SamplingTableRead)
    ReadSamplingTable(samplingTableFile);
  in.CloseSection();
  return true;
}


inline string DavidPAClass::SkipTo(ifstream &infile,string skipToString)
{
  string lineString;
  do{
    getline(infile,lineString);

  } while (lineString.find(skipToString,0)==string::npos && !infile.eof());

  if (infile.eof()){
    cerr<<"ERROR!!!  No "<<skipToString<<" found in Davids squarer file\n";
  }
  return lineString;

}


inline void DavidPAClass::ReadFORTRAN3Tensor(ifstream &infile,blitz::Array<double,3> &tempUkj)
{

   
  for (int counterTau=0;counterTau<tempUkj.extent(2);counterTau++){
    for (int counterUkj=0;counterUkj<tempUkj.extent(1);counterUkj++){
      for (int counterGridPoints=0;counterGridPoints<tempUkj.extent(0);counterGridPoints++){

	infile>>tempUkj(counterGridPoints,counterUkj,counterTau);
	//	tempUkj(counterGridPoints,counterUkj,counterTau)=0;
      }
    }
  }

}


inline bool IsDigit(char c)
{
  return (((c>='0') && (c<='9')) || (c == '-'));
}

inline bool IsNumberChar(char c)
{
  bool IsNumber = false;
  
  IsNumber = IsNumber || ((c>='0')&&(c<='9'));
  IsNumber = IsNumber || (c=='-');
  IsNumber = IsNumber || (c=='.');
  return IsNumber;
}



inline int GetNextInt(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsDigit(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  int num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}


inline string GetNextWord(string &s)
{
  int counter=0;
  while (s[counter] == ' '){
    counter++;
  }
  int startSpot=counter;
  while (s[counter] != ' '){
    counter++;
  }
  string toReturn=s.substr(startSpot,counter-startSpot);
  s=s.substr(counter,s.size()-counter);
  return toReturn;
}

inline double GetNextDouble(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsNumberChar(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  double num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}

inline void DavidPAClass::Print()
{ cerr<<"hi\n";}
            



#endif
