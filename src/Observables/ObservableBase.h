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

#ifndef OBSERVABLE_BASE_H
#define OBSERVABLE_BASE_H


#include "../PathDataClass.h"
#include "../Common.h"
#include "ObservableVar.h"
#include "../EventClass.h"
#include "../IO/IO.h"
#include <algorithm>
#include <numeric>
#include "ctime"
#include "sys/time.h"
#include <tr1/unordered_map>

using namespace IO;


/// This is the parent class for all observables.  It contains
/// a pointer to PathData.
class ObservableClass : public EventClass
{
 protected:

  template <typename Container> // we can make this generic for any container [1]
  struct container_hash {
    std::size_t operator()(Container const& c) const {
      double tot = 0;
      double sqtot = 0;
      for (int i = 0; i < c.size(); i++) {
        tot += c[i];
        sqtot += c[i]*c[i];
      }
      double value = sqtot/tot;
      return (size_t)*(const unsigned int *)&value;
    }
  };

  class CompareVecInt
  {
    public:
      bool operator() (const vector<int> &a, const vector<int> &b) {
        for (int i = 0; i<a.size(); i++)
          if (a[i] != b[i])
            return (a[i] > b[i]);
        return (a[0]>b[0]);
      }
  };

  //vector< vector<int> > PossPerms;
  //std::tr1::unordered_map<vector<int>,int,container_hash< vector<int> > > PossPerms;
  //std::tr1::unordered_map<vector<int>,int,container_hash< vector<int> > >::const_iterator PossPermsIterator;
  std::map<vector<int>,int,CompareVecInt> PossPerms;
  std::map<vector<int>,int,CompareVecInt>::const_iterator PossPermsIterator;
  Array<bool,1> CountedAlready;
  Array<int,1> TotalPerm;
  void SetupPermSectors(int n, int MaxNSectors=0);
  void GetPermInfo(vector<int> &ThisPerm, int &PermSector, int &PermNumber);

  bool TrackSign;



  struct CompareVectors
  {
    inline bool operator() (const vector<int> &a, const vector<int> &b) {
      //int amax = *max_element(a.begin(),a.end());
      //int bmax = *max_element(b.begin(),b.end());
      //if (amax == bmax)
      //  return b.size() < a.size();
      //else
      //  return amax < bmax;
      double atot = 0;
      double asqtot = 0;
      for (int i=0; i<a.size(); i++) {
        atot += a[i];
        asqtot += a[i]*a[i];
      }
      double btot = 0;
      double bsqtot = 0;
      for (int i=0; i<b.size(); i++) {
        btot += b[i];
        bsqtot += b[i]*b[i];
      }
      return (bsqtot/btot) > (asqtot/atot);

    }
  };

 public:
  /// The first time you write to an observable you have to do the
  /// write a little differently and you might need to write additional
  /// info like the description, etc.
  bool FirstTime;
  double SecondsInObservable;
  /// This stores how often this observable accumulates.  E.g. if
  /// Frequency is 3, Accumulate will actually accumulate every third
  /// time it is encountered the the algorithm.
  int Frequency;
  /// Same as Frequency, only getting track of time
  int TemporalFrequency;
  struct timeval starttime;
  struct timezone tzone;
  /// This a convenience function that allows one to specify a unit
  /// conversion if desired.  Set to 1.0 by default.
  double Prefactor;
  /// Note: This is not a reference.  If it were, it could change
  /// behind our backs
  string Description;
  /// Observe the state of the present path and add it to the
  /// running sum for averages.
  virtual void Accumulate() = 0;
  virtual void WriteBlock()=0;
  virtual void Read(IOSectionClass& IO);
  virtual void WriteInfo();

  /// This will just call Accumulate() every Frequency time it is called.
  void DoEvent();

  /// Calculate weight for sign tracking and/or importance sampling
  double CalcFullWeight();

  /// The constructor.  Sets PathData references and calls initialize.
  /// Note: the ioSection is passed by value, NOT by reference.  This
  /// is so we maintain our position in the file, even if the section
  /// variable we pass here changes.  If we pass the IO section to 
  /// any ObservableVar classes, we should do with our local IOSection
  /// variable, NOT the reference passed to derived classes.
  ObservableClass(PathDataClass &pathData,IOSectionClass &out) 
    : EventClass (pathData, out), FirstTime(true), Prefactor(1.0), TrackSign(false)
  {
  }

};




#endif
