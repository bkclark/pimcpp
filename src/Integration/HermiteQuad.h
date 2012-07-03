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

#ifndef HERMITE_QUAD_H
#define HERMITE_QUAD_H

#include <cmath>
#include <blitz/array.h>

using namespace blitz;

class Hermite16
{
public:
  static const int n=16;
  static const double u[16];
  static const double weight[16];
  double w[16];
  double x[16];
};

class Hermite20
{
public:
  static const int n=20;
  static const double u[20];
  static const double weight[20];
  double w[20];
  double x[20];
};

class Hermite30
{
public:
  static const int n=30;
  static const double u[30];
  static const double weight[30];
  double w[30];
  double x[30];
};



template <class RuleClass, class ScaledIntegrand> 
class HermiteQuadClass
{
 private:
  RuleClass Rule;
  double sigma;
 public:
  void SetSigma(double mysigma)
  {
    sigma = mysigma;
    for (int i=0; i<Rule.n; i++)
      Rule.x[i] = Rule.u[i]*M_SQRT2*sigma;
  }  
  
  inline double Integrate(ScaledIntegrand &Integrand)
  {
    double sum=0; 
    for (int i=0; i<Rule.n; i++)
      sum += Rule.w[i]*Integrand(Rule.x[i]);
    sum *= M_SQRT2*sigma;
    return (sum);
  }
  HermiteQuadClass()
  {
    for (int i=0; i<Rule.n; i++)
      Rule.w[i] = Rule.weight[i]*exp(-Rule.u[i]*Rule.u[i]);
  }
};



template <class RuleClass, class ScaledIntegrand> 
class Hermite3DQuadClass
{
 private:
  RuleClass Rule;
  double sigma;
  Array<double,2> Points;
 public:
  void SetSigma(double mysigma) 
  {
    int n = Rule.n;
    Points.resize(n*n*n,4);
    
    sigma = mysigma;
    for (int i=0; i<Rule.n; i++)
      Rule.x[i] = Rule.u[i]*M_SQRT2*sigma;
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
	for (int k=0; k<n; k++) {
	  int index = k+(j+i*n)*n;
	  Points(index,0) = Rule.x[i];
	  Points(index,1) = Rule.x[j];
	  Points(index,2) = Rule.x[k];
	  Points(index,3) = Rule.w[i]*Rule.w[j]*Rule.w[k];
	}
  }  
  
  inline double Integrate(ScaledIntegrand &Integrand)
  {
    int n = Rule.n;
    int N = n*n*n;
    double sum=0; 
    for (int i=0; i<N; i++)
      sum += Points(i,3)*Integrand(Points(i,0), Points(i,1), Points(i,2));
    double a = M_SQRT2*sigma;
    sum *= (a*a*a);
    return (sum);
  }

  Hermite3DQuadClass()
  {
    for (int i=0; i<Rule.n; i++)
      Rule.w[i] = Rule.weight[i]*exp(-Rule.u[i]*Rule.u[i]);
  }
};

#endif
