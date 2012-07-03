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

#include "RungeKutta.h"
#include "Integrate.h"

class SHO
{
public:
  double k;
  inline Vec2 operator()(double t, const Vec2 &x_and_dx)
  {
    Vec2 deriv;
    deriv[0] = x_and_dx[1];
    double x = x_and_dx[0];
    deriv[1] = -k*x;
    //    cerr << "deriv = [" << deriv[0] << ", " << deriv[1] << "]\n";
    return(deriv);
  }
};


void TestRK()
{
  SHO sho;
  sho.k = 1.0;
  LinearGrid t(0, 25.0, 1001);
  Array<Vec2,1> x(t.NumPoints);
  RungeKutta<SHO,Vec2> Integrator(sho);
  x(0)[0] = 0;
  x(0)[1] = 1.0;
  Integrator.Integrate(t, 0, t.NumPoints-1, x);
  for (int i=0; i<t.NumPoints-1; i++) 
    fprintf (stderr, "%1.12e %1.12e %1.12e\n", t(i), x(i)[0], x(i)[1]);
}

void RKSpeed()
{
  SHO sho;
  sho.k = 1.0;
  LinearGrid t(0, 25.0, 5001);
  Array<Vec2,1> x(t.NumPoints);
  RungeKutta<SHO,Vec2> Integrator(sho);
  x(0)[0] = 0;
  x(0)[1] = 1.0;
  clock_t start = clock();
  for (int i=0; i<10000; i++)
    Integrator.Integrate(t, 0, t.NumPoints-1, x);
  clock_t end = clock();
  double speed = 10000.0/(double)(1e-6*(end-start));
  cerr << "Runge Kutta: " << speed << " integrations per second.\n";
}

Vec2 ShoDeriv (double t, Vec2 u_and_du, void *shoPtr)
{
  SHO &sho = *((SHO *)shoPtr);
  return (sho(t, u_and_du));
}

void IntegrateSpeed()
{
  SHO sho;
  sho.k = 1.0;
  LinearGrid t(0, 25.0, 5001);
  Array<Vec2,1> x(t.NumPoints);
  x(0)[0] = 0;
  x(0)[1] = 1.0;

  clock_t start = clock();
  for (int i=0; i<10000; i++)
    IntegrateSecondOrder (t, 0, t.NumPoints-1, x,
			  ShoDeriv, &sho);
  clock_t end = clock();
  double speed = 10000.0/(double)(1e-6*(end-start));
  cerr << "Integrate: " << speed << " integrations per second.\n";
}





main()
{
  RKSpeed();
  IntegrateSpeed();
  //  TestRK();
}
