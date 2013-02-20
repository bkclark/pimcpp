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

#ifndef RANDOM_CLASS
#define RANDOM_CLASS
#include <sprng.h>
#include "../Communication/Communication.h"
#include <math.h>

class RandomClass
{
private:
  int NumClones;
  int CloneNumber;
  int *LocalStream, *CommonStream, *WorldStream;
  CommunicatorClass &MyComm;
  int NumCommon, NumLocal, NumWorld;
public:
  ///Produces a double that is unique to your processor
  inline double Local()
  { NumLocal++;  return sprng(LocalStream);  }

  ///Produces a double that is common among this clone
  inline double Common()
  { NumCommon++; return sprng(CommonStream); }
    
  ///Produces a double that is common to all processors
  inline double World()
  { NumWorld++;  return sprng(WorldStream);  }
  
  ///Produces an int that is unique to your processor between 0 and max
  inline int LocalInt(int max)
  { return (int)floor(Local()*max); }
  inline int CommonInt(int max)
  { return (int)floor(Common()*max); }
  inline int WorldInt(int max)
  { return (int)floor(World()*max); }
  
  /// normal random variate generator 
  /// Returns a random number distributed according to
  /// P(x) = 1/sqrt(2*pi * sigma) * exp (x^2 /(2*sigma^2))
  inline double LocalGaussian  (double sigma);
  inline double LocalGaussian2 (double sigma);
  inline double CommonGaussian (double sigma);
  inline double WorldGaussian  (double sigma);
  
  /// Produces a guassian random vector with radius that has STD sigma
  inline void LocalGaussianVec (double sigma, Vec1 &c);
  inline void CommonGaussianVec(double sigma, Vec1 &c);
  inline void WorldGaussianVec (double sigma, Vec1 &c);


  inline void LocalGaussianVec (double sigma, Vec2 &c);
  inline void CommonGaussianVec(double sigma, Vec2 &c);
  inline void WorldGaussianVec (double sigma, Vec2 &c);

  inline void LocalGaussianVec (double sigma, Vec3 &c);
  inline void CommonGaussianVec(double sigma, Vec3 &c);
  inline void WorldGaussianVec (double sigma, Vec3 &c);
      
  void Init()
  {
    int seed = make_sprng_seed();
    //cerr<<"My seed is "<<seed<<endl;
    Init (seed);
  }
  
  int make_new_seed()
  {
    
    time_t tp;
    struct tm *temp;
    unsigned int temp2, temp3;
    static unsigned int temp4 = 0xe0e1;
    
    time(&tp);
    temp = localtime(&tp);
    
    temp2 = (temp->tm_sec<<26)+(temp->tm_min<<20)+(temp->tm_hour<<15)+
      (temp->tm_mday<<10)+(temp->tm_mon<<6);
    temp3 = (temp->tm_year<<13)+(temp->tm_wday<<10)+(temp->tm_yday<<1)+
      temp->tm_isdst;
    temp2 ^= clock()^temp3;
    
    temp4 = (temp4*0xeeee)%0xffff;
    temp2 ^= temp4<<16;
    temp4 = (temp4*0xaeee)%0xffff;
    temp2 ^= temp4;
    
    temp2 &= 0x7fffffff;
    
    return temp2;
  }
  
  int InitWithRandomSeed(int numClones=1)
  {
    int seed;
    if (MyComm.MyProc()==0){
      seed = make_new_seed();
      seed = 123456;
      cout<<"THE RANDOM SEED IS "<<seed<<endl;
    }
    MyComm.Broadcast(0,seed);
    Init (seed, numClones);
    return seed;
  }
  
  
  void Init(int seed, int numClones=1, bool SameSeed=false)
  {
    NumClones=numClones;
    int myProc=MyComm.MyProc();
    int numProcs=MyComm.NumProcs();
    int procsPerClone=numProcs/NumClones;
    CloneNumber = myProc/procsPerClone;
    int commID  = numProcs+CloneNumber;
    int worldID = numProcs+NumClones;
    if(!SameSeed) {
      ///LocalStream is a stream of random numbers unique to the node
      LocalStream = 
        init_sprng(SPRNG_DEFAULT,myProc,numProcs+NumClones+1,seed,SPRNG_DEFAULT);
    } else {
      ///unless you explicitly sync the streams
      LocalStream = 
        init_sprng(SPRNG_DEFAULT,0,numProcs+NumClones+1,seed,SPRNG_DEFAULT);
    }
    /// CommonStream is a stream of random numbers shared between all
    /// processors that are on the same clone 
    CommonStream = 
      init_sprng(SPRNG_DEFAULT,commID,numProcs+NumClones+1,seed,SPRNG_DEFAULT);
    /// World Stream is a stream of random numbers shared between all
    /// processors in MPICommWorld
    WorldStream =  init_sprng(SPRNG_DEFAULT,worldID,numProcs+NumClones+1,
			      seed,SPRNG_DEFAULT);
  }

  RandomClass(CommunicatorClass &comm) : 
    MyComm(comm), NumCommon(0), NumLocal(0)
  {
    
  }
};

inline double RandomClass::LocalGaussian(double sigma)
{                                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    do {
      V1 = 2.0 * Local() - 1.0;
      V2 = 2.0 * Local() - 1.0;
      S = V1 * V1 + V2 * V2;
    } while ( S >= 1.0 );
    
    temp = sqrt( (-2.0 * log((double) S ) ) / S );
    X1 = V1 * temp;
    X2 = V2 * temp;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}


inline double RandomClass::LocalGaussian2(double sigma)
{                                      
  double X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    double v1 = Local();
    double v2 = Local();
    double mag = sqrt(-2.0*log(v1));
    double s, c;
#ifdef HAVE_SINCOS
    sincos(2.0*M_PI*v2, &s, &c);
#else
    s = sin(2.0*M_PI*v2);
    c = cos(2.0*M_PI*v2);
#endif
    X1 = mag * c;
    X2 = mag * s;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}


inline double 
RandomClass::CommonGaussian(double sigma)
{                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    do {
      V1 = 2.0 * Common() - 1.0;
      V2 = 2.0 * Common() - 1.0;
      S = V1 * V1 + V2 * V2;
    } while ( S >= 1.0 );
    
    temp = sqrt( (-2.0 * log((double) S ) ) / S );
    X1 = V1 * temp;
    X2 = V2 * temp;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}

inline double 
RandomClass::WorldGaussian(double sigma)
{                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    do {
      V1 = 2.0 * World() - 1.0;
      V2 = 2.0 * World() - 1.0;
      S = V1 * V1 + V2 * V2;
    } while ( S >= 1.0 );
    
    temp = sqrt( (-2.0 * log((double) S ) ) / S );
    X1 = V1 * temp;
    X2 = V2 * temp;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}

inline void 
RandomClass::LocalGaussianVec(double sigma, Vec3 &c)
{ 
  c(0)=LocalGaussian(sigma); 
  c(1)=LocalGaussian(sigma); 
  c(2)=LocalGaussian(sigma); 
}

inline void 
RandomClass::CommonGaussianVec(double sigma, Vec3 &c)
{ 
  c(0)=CommonGaussian(sigma); 
  c(1)=CommonGaussian(sigma); 
  c(2)=CommonGaussian(sigma); 
}

inline void 
RandomClass::WorldGaussianVec(double sigma, Vec3 &c)
{ 
  c(0)=WorldGaussian(sigma); 
  c(1)=WorldGaussian(sigma); 
  c(2)=WorldGaussian(sigma); 
}

/// Produces a guassian random vector with radius that has STD sigma 
inline void 
RandomClass::LocalGaussianVec(double sigma,Vec2 &c)
{ 
  c(0)=LocalGaussian(sigma); 
  c(1)=LocalGaussian(sigma); 
}
 
inline void 
RandomClass::CommonGaussianVec(double sigma,Vec2 &c)
{ 
  c(0)=CommonGaussian(sigma); 
  c(1)=CommonGaussian(sigma); 
}

inline void 
RandomClass::WorldGaussianVec(double sigma,Vec1 &c)
{ 
  c(0)=WorldGaussian(sigma); 
}

/// Produces a guassian random vector with radius that has STD sigma 
inline void 
RandomClass::LocalGaussianVec(double sigma,Vec1 &c)
{ 
  c(0)=LocalGaussian(sigma); 
}
 
inline void 
RandomClass::CommonGaussianVec(double sigma,Vec1 &c)
{ 
  c(0)=CommonGaussian(sigma); 
}

inline void 
RandomClass::WorldGaussianVec(double sigma,Vec2 &c)
{ 
  c(0)=WorldGaussian(sigma); 
  c(1)=WorldGaussian(sigma); 
}



#endif
