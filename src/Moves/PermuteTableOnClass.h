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

//#include "../Common/Blitz.h"
#include "../PathDataClass.h"
#ifndef PERMUTE_TABLE_ON_CLASS_H
#define PERMUTE_TABLE_ON_CLASS_H

// class Hclass
// {
//  public:
//   int j;
//   double Htilde_ij;
// };

/// Note:  all indices in this file are with respect to the present
/// species only.  Thus, when accessing the actual path, we must add
/// the PathData.Species(SpeciesNum).FirstPtcl to our indices.

// class CycleClass
// {
//  public:
//   int Length;
//   TinyVector<int,4> CycleRep;
//   double P, C;
//   ///Takes the number of particles you have and returns the
//   ///representation such that p(j) is the particle that the j'th
//   ///particle permutes onto. 
//   //  void CanonicalPermRep(Array<int,1> myArray);
//   void Apply (PathClass &path, int firstPtcl, int timeSlice);

// };
// inline bool operator==(const CycleClass &a,const CycleClass &b)
// {
//   if (a.Length!=b.Length){
//     return false;
//   }
//   for (int counter=0;counter<a.Length;counter++){
//     if (a.CycleRep[counter]!=b.CycleRep[counter])
//       return false;
//   }
//   return true;
// }


struct ParticleProbStruct
{
  int ptcl;
  double prob;


};


class PermuteTableOnClass
{
 private:
  int TableSize;
  int SpeciesNum;
  int Slice1, Slice2;
  inline void AddEntry(const CycleClass &cycle);
  int ExcludeParticle;
  PathDataClass &PathData;
  void ConstructHTable();
  double calcij(int i, int j);
  double Norm, NormInv;
  inline int FindEntry(double xi);  
  inline int FindEntrySlow(double xi);  
public:
  ///Designed to allow it to give only odd or only even permutations.  Specifically used for a coupling move.
  bool OnlyOdd;
  bool OnlyEven;
  int NumEntries;
  CycleClass CurrentCycle;
  Array<list<ParticleProbStruct>,1> HTable;
  //  Array<double,2> HTable;
  Array<CycleClass,1> CycleTable;

  //  void PermuteHTable();

  /// The enhancement factor for each cycle length. That is we
  /// multiply the probability of a 3-cycle by Gamma[2]  
  TinyVector<double,4> Gamma; 
  
  // Smallest value of exp(-s^2/(4*lambda*beta) which we will include
  // in Htable.
  double epsilon;

  void ConstructCycleTable(int speciesNum,int slice1,int slice2);
  void ConstructCycleTable(int speciesNum,int slice1,int slice2,
			   int excludeParticle);
  void CanonicalPermRep(Array<int,1> P);
  double AttemptPermutation();
  double CalcReverseProb(const PermuteTableOnClass &forwardTable);
  Array<int,1> CurrentParticles();



  void Read(IOSectionClass &inSection);

  PermuteTableOnClass(PathDataClass &myPathData) : PathData(myPathData)
  {
    NumEntries=0;
    TableSize = 1000;
    CycleTable.resize(TableSize);
    OnlyOdd=false;
    OnlyEven=false;
    ExcludeParticle=-1;
  }

};

inline void PermuteTableOnClass::AddEntry(const CycleClass &cycle)
{
  if (NumEntries >= (TableSize-1)) {
    TableSize *= 2;
    cerr<<"My new tablesize is "<<TableSize<<endl;
    CycleTable.resizeAndPreserve(TableSize);
  }

  CycleTable(NumEntries) = cycle;
  NumEntries++;
}


inline int PermuteTableOnClass::FindEntrySlow(double xi)
{
  xi *= Norm; 
  int num=0;
  while (CycleTable(num).C<xi){
    num++;
  }
  return num;
  
}  

//Pass a random number between 0 and 1 to this function and it
//returns the index of the permutation with the appropriate
//probability. 
inline int PermuteTableOnClass::FindEntry(double xi) 
{
  // Do a binary search
  //  int toCheck=FindEntrySlow(xi);
  xi *= Norm; 
  int hi = NumEntries-1;
  int lo = 0;

  if (xi < CycleTable(0).C){
//     if (0!=toCheck){
//       cerr<<"ERROR! ERROR! We are DUMB!"<<endl;
//     }
    return (0);
  }


  while (hi-lo>1){
    int attempt = (hi+lo)>>1;
    if (xi<CycleTable(attempt).C)
      hi=attempt;
    else 
      lo=attempt;
  }
//   if (hi!=toCheck)
//       cerr<<"ERROR! ERROR! We are DUMB!"<<endl;
  return hi;
}


//   while (attempt != lo) {
//     attempt = (hi+lo)/2;
//     if (CycleTable(attempt).C > xi)
//       hi = attempt;
//     else
//       lo = attempt;
//   }
//   if (hi!=toCheck){
//     cerr<<"ERROR! ERROR! We are DUMB!"<<endl;
//   }
//   return (hi);
  
//   int num=0;
//   while (CycleTable(num).C<xi){
//     num++;
//   }
//   return num;



#endif
