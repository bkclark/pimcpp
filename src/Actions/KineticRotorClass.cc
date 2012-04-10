////////////////////////////////////////////////////////////
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

#include "../PathDataClass.h"
#include "KineticRotorClass.h"
#include "../Moves/MoveUtils.h"

KineticRotorClass::KineticRotorClass(PathDataClass &pathData ) : 
  RotorActionBaseClass (pathData)
{
}

///This has to be called after pathdata knows how many
///particles it has
void KineticRotorClass::Read(IOSectionClass& in)
{
  // read in moments of inertia
  // these are the defaults for water in Hartree(?)
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  in.ReadVar("Ia",Ia);
  in.ReadVar("Ib",Ib);
  in.ReadVar("Ic",Ic);
  A = 1.0/(2*Ia);//1.1411e-4;
  B = 1.0/(2*Ib);//5.7207e-5;
  C = 1.0/(2*Ic);//3.8105e-5;

  //int setJmax;
  //assert(in.ReadVar("SetJmax",setJmax));
  //int setPrec = 100;
  //assert(in.ReadVar("SetPrecision",setPrec));
  //rho = new RotorRhoClass(A, B, C, setJmax, setPrec); 
  //rho->SetTau(PathData.Path.tau);
  ////rho->Init(setJmax, 0., M_PI/2, 1000);
  //rho->Init(setJmax, 0., M_PI, 1999);
  //cerr << "WARNING: FOR TESTING I AM doing explicit comutation of rho: THIS WILL BE SLOW" << endl;

  ReadGridPoints(in);
  ReadEnergyGridPoints(in);
  doESq = false;
  in.ReadVar("ReadEnergySquared",doESq);
  if(doESq) {
    string esqfilename;
    assert(in.ReadVar("EnergySquaredOutput",esqfilename));
    esqOut.open(esqfilename.c_str());
  }
  cerr << "LEAVING KINETIC ROTOR READ" << endl;

  // debugging init
  //err.open("Rotor.err.dat");
  //tolrho.open("Rotor.rho.tolerance.dat");
  //overX.open("Rotor.overflow.exact.dat");
  //overS.open("Rotor.overflow.spline.dat");
  thist.resize(300);
  phist.resize(300);
  chist.resize(300);
  for(int z=0; z<thist.size(); z++) {
    thist(z) = 0;
    phist(z) = 0;
    chist(z) = 0;
  }
  cerr << "arrays init " << thist.size() << " " << phist.size() << " " << chist.size() << endl;
  DT = M_PI/300;
  DP = 2*M_PI/300;
  times=0;
  maxRhoS = 0.;
  maxRhoX = 0.;
  TotalZ = 0.; TotalE = 0.; TotalESq = 0.;
  THETA = 0; PHI = 0; CHI = 0;
}

void KineticRotorClass::ReadGridPoints(IOSectionClass& in)
{
  // read in num of grid points
  int Ntheta, Nphi, Nchi;
  assert(in.ReadVar("NumThetaPoints",Ntheta));
  assert(in.ReadVar("NumPhiPoints",Nphi));
  bool oldFormat = false;
  in.ReadVar("OldFormat",oldFormat);
  //assert(in.ReadVar("NumChiPoints",Nchi));
  //double last = M_PI * (1. - 1./(Ntheta+1));
  //ThetaGrid = new LinearGrid(0., last, Ntheta);

  double lastTheta = M_PI * (1. - 1./Ntheta);
  double lastPhi = 4*M_PI * (1. - 1./Nphi);
  if(oldFormat) {
    lastTheta = M_PI;
    lastPhi = 4 * M_PI;
  }
  //PhiGrid = new LinearGrid(0., last, Nphi);
  //last = 2*M_PI * (1. - 1./(Nchi+1));
  //ChiGrid = new LinearGrid(0., last, Nchi);

  // hack testing spline init
  //ThetaGrid = new LinearGrid(0., M_PI, Ntheta);
  ThetaGrid = new LinearGrid(0., lastTheta, Ntheta);
  //PhiGrid = new LinearGrid(0., 2*M_PI, Nphi);
  //PhiGrid = new LinearGrid(0., 4*M_PI, Nphi);
  PhiGrid = new LinearGrid(0., lastPhi, Nphi);
  //ChiGrid = new LinearGrid(0., 2*M_PI, Nchi);

  // read in spline gridpoints
  string filename;
  assert(in.ReadVar("File",filename));
  //Array<double,3> values;
  Array<double,2> values;
  values.resize(Ntheta, Nphi);
  double vol = Ntheta * Nphi;// * Nchi;
  int count = 0;
  ifstream infile(filename.c_str());
  cerr << "Reading in grid points from" << filename << endl;
  string line;
  cerr << "Getting Jmax: ";
  infile >> line;
  cerr << line << endl;
  cerr << "Getting tau: ";
  infile >> line;
  cerr << line << endl;
  //while(infile) {}
  int countT=0;
  for(int nt=0; nt<Ntheta; nt++) {
    for(int np=0; np<Nphi; np++) {
      //for(int nc=0; nc<Nchi; nc++) {
      double theta, phi, chi, logRho;
      infile >> line;
      theta = atof(line.c_str());
      infile >> line;
      phi = atof(line.c_str());
      //infile >> line;
      //chi = atof(line.c_str());
      infile >> line;
      double checkRho = atof(line.c_str());
      infile >> line;
      infile >> line;

      // checking readin of angles
      double checkT, checkP;//, checkC;
      checkT = (*ThetaGrid)(nt);
      checkP = (*PhiGrid)(np);
      //checkC = (*ChiGrid)(nc);
      //cerr << "Comparing angles at " << nt << " " << np << endl;//" " << nc << endl;
      //cerr << "theta " << theta << " " << checkT << endl;
      //cerr << "phi " << phi << " " << checkP << endl;
      //cerr << "chi " << chi << " " << checkC << endl;
      assert(abs(checkT - theta) < 1e-5);
      assert(abs(checkP - phi) < 1e-4);
      //assert(abs(checkC - chi) < 1e-5);

      logRho = atof(line.c_str());
      //values(nt, np, nc) = logRho;
      values(nt, np) = logRho;
      if(checkP == 0) {//) && checkC == 0) {
        countT++;
        //cout << countT << " " << (*ThetaGrid)(nt) << " " << (*PhiGrid)(np) << " " << values(nt,np) << " " << exp(values(nt,np)) << endl;
      }
      count ++;
      //}
    }
  }

  cerr << "Read in " << count << " line; expected " << vol << endl;
  cerr << "Array extents " << values.extent(0) << " " << values.extent(1) << endl;
  cerr << "Grid sizes " << ThetaGrid->NumPoints << " " << PhiGrid->NumPoints << endl;
  spline.Init(ThetaGrid, PhiGrid, values);
  cerr << "Spline initialized" << endl;
}

void KineticRotorClass::ReadEnergyGridPoints(IOSectionClass& in)
{
  // read in num of grid points
  int Ntheta, Nphi, Nchi;
  assert(in.ReadVar("NumEnergyThetaPoints",Ntheta));
  assert(in.ReadVar("NumEnergyPhiPoints",Nphi));
  //assert(in.ReadVar("NumEnergyChiPoints",Nchi));
  bool oldFormat = false;
  in.ReadVar("EnergyOldFormat",oldFormat);
  double lastTheta = M_PI * (1. - 1./Ntheta);
  double lastPhi = 4*M_PI * (1. - 1./Nphi);
  if(oldFormat) {
    lastTheta = M_PI;
    lastPhi = 4 * M_PI;
  }
  ThetaEnergyGrid = new LinearGrid(0., lastTheta, Ntheta);
  PhiEnergyGrid = new LinearGrid(0., lastPhi, Nphi);
  //ThetaEnergyGrid = new LinearGrid(0., M_PI, Ntheta);
  //PhiEnergyGrid = new LinearGrid(0., 4*M_PI, Nphi);
  //ChiEnergyGrid = new LinearGrid(0., 2*M_PI, Nchi);

  // read in spline gridpoints
  string filename;
  assert(in.ReadVar("EnergyFile",filename));
  Array<double,2> values;
  values.resize(Ntheta, Nphi);
  double vol = Ntheta * Nphi;
  int count = 0;
  ifstream infile(filename.c_str());
  cerr << "Reading in Energy grid points from" << filename << endl;
  string line;
  cerr << "Getting Jmax: ";
  infile >> line;
  cerr << line << endl;
  cerr << "Getting tau: ";
  infile >> line;
  cerr << line << endl;
  for(int nt=0; nt<Ntheta; nt++) {
    for(int np=0; np<Nphi; np++) {
      //for(int nc=0; nc<Nchi; nc++) {
      double theta, phi, chi, E;
      infile >> line;
      theta = atof(line.c_str());
      infile >> line;
      phi = atof(line.c_str());
      infile >> line;
      chi = atof(line.c_str());
      infile >> line;

      // checking readin of angles
      double checkT, checkP, checkC;
      checkT = (*ThetaEnergyGrid)(nt);
      checkP = (*PhiEnergyGrid)(np);
      //checkC = (*ChiEnergyGrid)(nc);
      //cerr << "Comparing angles" << endl;
      //cerr << "theta " << theta << " " << checkT << endl;
      //cerr << "phi " << phi << " " << checkP << endl;
      //cerr << "chi " << chi << " " << checkC << endl;
      assert(abs(checkT - theta) < 1e-5);
      assert(abs(checkP - phi) < 1e-4);
      //assert(abs(checkC - chi) < 1e-5);

      E = atof(line.c_str());
      //values(nt, np, nc) = E;
      values(nt, np) = E;
      //cerr << (*ThetaEnergyGrid)(nt) << " " << (*PhiEnergyGrid)(np) << " " << (*ChiEnergyGrid)(nc) << values(nt,np,nc) << endl;
      count ++;
      //}
    }
  }

  cerr << "Read in " << count << " line; expected " << vol << endl;
  //EnergySpline.Init(ThetaEnergyGrid, PhiEnergyGrid, ChiEnergyGrid, values);
  EnergySpline.Init(ThetaEnergyGrid, PhiEnergyGrid, values);
  cerr << "Energy Spline initialized" << endl;
}

double 
KineticRotorClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles, int level)
{
  //cerr << "CONVENTIONAL KineticRotor slices " << slice1 << " " << slice2 << endl;
  //cerr << "KineticRotor testing..." << endl;
  //int numP = 1000;
  ////double dt = M_PI/numP;
  //double dt = 0.00314;
  ////double t0 = dt*numP;
  //double t0 = 0.;
  //for(int t=0; t<numP; t++) {
  //  double theta, phi, chi;
  //  theta = dt*t + t0;
  //  phi = 0.;
  //  chi = 0.;
  //  //for(int p=0; p<numP; p++) {
  //  //  phi = 2*(dt*p + t0);
  //  //  for(int c=0; c<numP; c++) {
  //  //    chi = 2*(dt*c + t0);
  //      double KE = EnergySpline(theta, phi + chi);
  //      double logR = spline(theta, phi + chi);
  //      cout << theta << " " << phi << " " << chi << " " << logR << " " << exp(logR) << " " << KE << endl;
  //  //  }
  //  //}
  //}
  //exit(1);
  //////////////////////////////////

  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int numChangedMol = molRoster.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    //double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      //double tInv, pInv, cInv;
      //getAngles(slice+skip, mol, slice, mol, tInv, pInv, cInv);
      //cerr << "ANGLES " << phi << " " << theta << " " << chi << endl;
      //cerr << "       " << pInv << " " << tInv << " " << cInv << endl;

      // explicit computation:
      //cerr << "Calling CalcRho with angles " << phi << " " << theta << " " << chi << " between slices " << slice << " " << slice+skip << endl;
      //double r = rho->CalcRhoAllJ(phi, theta, chi);
      //////cerr << endl;
      ////// need to check for negative rho; use noisy correction
      //if(r<0) {
      //  //err << "Got NEGATIVE exact rho " << r << " " << phi << " " << theta << " " << chi << endl;
      //  //double a = -1 * PathData.Path.Random.Local() * r;
      //  //cout << " generated noise " << a << endl;
      //  //r = a;
      //  r = abs(r);
      //}
      //double XlogRho = log(r);
      //////cerr << "Computed exact rho " << r << " and log " << XlogRho;// << endl;

      //cerr << "Accessed spline at indices " << theta << " " << phi << " " << chi << " value ";
      // 2D spline: rho is a function of **sum** of phi and chi
      double logRho = spline(theta, phi + chi);
      //cout << XlogRho << " " << logRho << " " << logRho - XlogRho << endl;
      //double logRho = spline(theta, phi, chi);
      //cerr << "Spline logRho " << logRho << endl; 
      //if(abs(XlogRho-logRho) > 1e-2)
      //  tolrho << theta << " " << phi << " " << chi << " " << XlogRho << " " << logRho << endl;
      TotalK -= logRho;
      //TotalK += XlogRho;
      //if(theta > 1) {
      //  cout << slice << " " << theta << " " << chi+phi << " " << logRho << " " << exp(logRho) << endl;
      //}
      //cout << slice << " " << theta << " " << chi+phi << " " << logRho << " || ";
    }
  }
  //cout << "returning Action " << TotalK << endl << endl;
  return (TotalK);
}

            // cang checkangles
            // get coords; check getAngles
            //dVec u(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
            //dVec ref=cross(u,(0.,0.,1.));
            //double cchi, schi;
            //cchi = cos(chi); schi=sin(chi);
            //dVec perp = cross(u,ref);
            //dVec p = Renormalize(ref,cchi) + Renormalize(perp, schi);
            ////double x = 2 * (PathData.Path.Random.Local() - 0.5);
            ////double y = 2 * (sqrt(1-x*x)) * (PathData.Path.Random.Local() - 0.5);
            ////double z = (-x*u(0) -y*u(1))/u(2);
            ////chi = GetAngle(ref, p);
            //err << theta << " " << phi << " " << chi << " ";// << u << " " << p << endl;
            ////err << p(0) << " " << p(1) << " " << p(2) << " " << Mag(p) << " " << chi << endl;
            //u = Renormalize(u, 0.58588);
            //p = Renormalize(p, 0.75695);
            //dVec p2a = u + p;
            //dVec p2b = u - p;
            //dVec p1a(-0.75695, 0., 0.58588);
            //dVec p1b(0.75695, 0., 0.58588);

            //// FROM GETANGLES
            //// compute bisector, norm, and in-plane unit vectors
            //dVec b1 = GetBisector(p1a, p1b);
            //dVec n1 = Normalize(cross(p1a, p1b));
            //dVec u1 = Normalize(cross(b1, n1));
            //dVec b2 = GetBisector(p2a, p2b);
            //dVec n2 = Normalize(cross(p2a, p2b));
            //dVec u2 = Normalize(cross(b2, n2));

            //// theta is angle between in-plane vectors u
            //// (rotation about norm to molecule plane)
            //if(isEqual(u1, u2)) {
            //  theta = 0.;
            //  phi = GetAngle(b2, b1);
            //  chi = 0.;
            //}
            //else if(isEqualOpp(u1, u2)) {
            //  theta = M_PI;
            //  dVec newB = -1 * b1;
            //  phi = GetAngle(b2, newB);
            //  chi = 0.;
            //} else {
            //  theta = GetAngle(u2, u1);
            //  // axis of rotation for theta is cross of in-plane vectors u
            //  dVec A;
            //  A = Normalize(cross(u1,u2));

            //  // phi is angle between norm n1 and A
            //  phi = GetAngle(n1, A);
            //
            //  // chi is angle between norm A and n2
            //  chi = GetAngle(A, n2);
            //}
            //err << theta << " " << phi << " " << chi << endl;

double KineticRotorClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  //cerr << "COORDS SLICE0 SLICE1 SLICE2" << endl;
  //cerr << PathData.Path(0,0) << PathData.Path(1,0) << PathData.Path(2,0) << endl;
  //cerr << PathData.Path(0,1) << PathData.Path(1,1) << PathData.Path(2,1) << endl;
  //cerr << PathData.Path(0,2) << PathData.Path(1,2) << PathData.Path(2,2) << endl;
  //
  int numS = 0;
  int M = PathData.Path.NumTimeSlices()-1;
  //cerr << "Calculating energy over slices " << slice1 << " " << slice2 << endl;
  double TotalK = 0.0;
  double TotalXK = 0.0;
  double TotalESq = 0.0;
  double Z = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double theta, phi, chi;
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      numS++;
      double theta, phi, chi;
      //cerr << "getting angles...";
      getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      //cerr << "got " << theta << " " << phi+chi << endl;

      //double XKE = rho->CalcEnergyAllJ(phi, theta, chi);
      //TotalZ += 1;
      //TotalE += XKE;
      //cout << TotalE/TotalZ << " " << XKE << endl;
      //cerr << slice << " " << slice+skip << " " << theta << " " << phi << " " << chi << endl;
      //int count=0;
      //for(int t=0; t<=200; t++) {
      //  for(int p=0; p<=500; p++) {
      //    for(int c=0; c<=500; c++) {
      //      //double phi = p * 2 * M_PI/50;
      //      //double chi = c * 2 * M_PI/50;
      //      //double theta = t * M_PI/200;
      //theta = M_PI * PathData.Path.Random.Local();
      //phi = 2 * M_PI * PathData.Path.Random.Local();
      //chi = 2 * M_PI * PathData.Path.Random.Local();
      //double w = sin(theta);// * M_PI/200 * M_PI/25 * M_PI/25;
      //      double w = 1.;

      //      double Z, E, ESq;
      //      //rho->Quadrature(phi + chi, theta, Z, E, ESq);
      //      //cerr << theta << " " << phi+chi << " " << Z << " " << E << " " << ESq << endl;
      //      Z = exp(spline(theta, phi+chi));
      //      E = EnergySpline(theta, phi+chi);
      //      ESq = 0.;
      //      TotalZ += w*Z*Z;
      //      TotalE += w*Z*Z*E;
      //      TotalESq += w*Z*ESq;
      //      //count++;
      // test angular distribution
      times++;
      int index = int(theta/DT);
      //thist(index) += w*Z*Z;
      thist(index) += 1;
      index = int(phi/DP);
      //phist(index) += w*Z*Z;
      phist(index) += 1;
      index = int(chi/DP);
      //chist(index) += w*Z*Z;
      chist(index) += 1;
      //      if(count%500000 == 0)
      //        cerr << count << " " << TotalE/TotalZ << " " << TotalZ << endl;
      //    }
      //  }
      //}
       

      // write out angular distributions
      //if(times%5000000==0) {
      //  tout.open("Rotor.theta.dat");
      //  pout.open("Rotor.phi.dat");
      //  chiout.open("Rotor.chi.dat");
      //  for(int n=0; n<300; n++) {
      //    //tout << DT*(n+1) << " " << thist(n)/TotalZ << " " << double(thist(n))/times*sin(DT*(n+1)) << endl;
      //    //pout << DP*(n+1) << " " << phist(n)/TotalZ << endl;
      //    //chiout << DP*(n+1) << " " << chist(n)/TotalZ << endl;
      //    tout << DT*(n+1) << " " << thist(n)/times << " " << double(thist(n))/times*sin(DT*(n+1)) << endl;
      //    pout << DP*(n+1) << " " << phist(n)/times << endl;
      //    chiout << DP*(n+1) << " " << chist(n)/times << endl;
      //  }
      //  tout.close();
      //  pout.close();
      //  chiout.close();
      //  //exit(1);
      //}
      
       
      //cout << theta << " " << phi+chi << " " << TotalE/TotalZ << " " << E/Z << endl;
      //exit(1);
        
        
      //typo error: negatiave sign in tabulated energy, corrected here!!
      //double KE = -1 * EnergySpline(theta, phi + chi);
      //cerr << "accessing spline " << theta << " " << phi+chi << endl;
      double z = exp(spline(theta, phi + chi));
      //cerr << "  " << z << endl;
      //cerr << "accessing energySpline " << theta << " " << phi+chi << endl;
      double KE = EnergySpline(theta, phi + chi);
      //cerr << "  " << KE << endl;
      //TotalE += KE*z*z;
      //TotalK += KE*z*z;
      ////TotalXK += XKE;
      //TotalZ += z*z;
      TotalE += KE;
      //TotalXK += XKE;
      TotalZ += 1;
      if(z < 1e-5 && KE>z)
        TotalK += 0.;
        //cerr << slice << " " << mol << " angles " << phi << " " << theta << " " << chi << " KE " << KE << " " << z << " ratio " << KE/z << " accum " << TotalK << endl;
      else
        TotalK += KE/z;
      if(doESq) {
        TotalESq += KE/(z*z);
      }
      //cerr << slice << " " << KE << " " << TotalK << endl;
      //if(abs(XKE) > maxRhoX) {
      //  overX << theta << " " << phi << " " << chi << " " << slice << " " << XKE << endl;
      //  maxRhoX = abs(XKE);
      //}
      //if(abs(KE) > maxRhoS) {
      //  overS << theta << " " << phi << " " << chi << " " << slice << " " << KE << " " << KE/z << endl;
      //  maxRhoS = abs(KE);
      //}
    }
  }
  if(doESq)
    esqOut << TotalESq << " " << TotalK << endl;
  //cerr << "returning " << TotalK/M << " slices " << numS << " " << " divide by " << M << endl;
  //cerr << "==============  Z  KE/Z  XKE/Z " << endl;
  //cout << TotalE/TotalZ << " " << 2*TotalE/TotalZ << " " << TotalZ << " " << endl;
  //cout << TotalK << " " << TotalK/M << " " << TotalE/TotalZ << " " << TotalE/(M*TotalZ) << endl;//" " << theta << " " << phi+chi << endl;
  return (TotalK);
  //return (TotalE/TotalZ);
  //cerr << "leaving dBeta" << endl;
}


string
KineticRotorClass::GetName()
{
  return "KineticRotor";
}

const double TOL = 1e-5;

bool isEqual(dVec u, dVec v)
{
  bool eq = false;
  if (abs(u(0) - v(0))< TOL) {
    if (abs(u(1) - v(1)) < TOL) {
      if (abs(u(2) - v(2)) < TOL) {
        eq = true;
      }
    }
  }
  return eq;
}

bool isEqualOpp(dVec u, dVec v)
{
  bool eq = false;
  if (abs(u(0) + v(0)) < TOL) {
    if (abs(u(1) + v(1)) < TOL) {
      if (abs(u(2) + v(2)) < TOL) {
        eq = true;
      }
    }
  }
  return eq;
}

// untested trial version
//void RotorActionBaseClass::getAngles(int slice1, int mol1, int slice2, int mol2, double& theta, double& phi, double& chi)
//{
//  // get molecule coordinates WRT molecule center
//  Array<int,1> mol1Members, mol2Members;
//  PathData.Mol.MembersOf(mol1Members, mol1);
//  PathData.Mol.MembersOf(mol2Members, mol2);
//  dVec p1a = PathData.Path(slice1, mol1Members(1)) - PathData.Path(slice1, mol1Members(0));
//  dVec p1b = PathData.Path(slice1, mol1Members(2)) - PathData.Path(slice1, mol1Members(0));
//  dVec p2a = PathData.Path(slice2, mol2Members(1)) - PathData.Path(slice2, mol2Members(0));
//  dVec p2b = PathData.Path(slice2, mol2Members(2)) - PathData.Path(slice2, mol2Members(0));
//
//  // compute bisector, norm, and in-plane unit vectors
//  dVec b1 = GetBisector(p1a, p1b);
//  dVec n1 = Normalize(cross(p1a, p1b));
//  dVec u1 = Normalize(cross(b1, n1));
//  dVec b2 = GetBisector(p2a, p2b);
//  dVec n2 = Normalize(cross(p2a, p2b));
//  dVec u2 = Normalize(cross(b2, n2));
//
//  double theta1 = acos(b1(2));
//  double theta2 = acos(b2(2));
//  double phi1 = acos(b1(0)/sin(theta1));
//  double phi2 = acos(b2(0)/sin(theta2));
//  dVec perp1 = p1a - Renormalize(b1,0.75695);
//  dVec ref1 = Normalize(cross(b1,(0.,0.,1.)));
//  double chi1 = GetAngle(ref1, perp1);
//  dVec perp2 = p2a - Renormalize(b2,0.75695);
//  dVec ref2 = Normalize(cross(b2,(0.,0.,1.)));
//  double chi2 = GetAngle(ref2, perp2);
//
//  theta = theta2 - theta1;
//  phi = phi2 - phi1;
//  chi = chi2 - chi1;
//  while(theta<0)
//    theta += M_PI;
//  while(phi<0)
//    phi += 2*M_PI;
//  while(chi<0)
//    chi += 2*M_PI;
//}

// map cartesian coordinates of two WATER molecules
// to an angular displacement expressed in terms of Euler angles
// Note this is specialized for a 3D rotor such as WATER!!
void RotorActionBaseClass::getAngles(int slice1, int mol1, int slice2, int mol2, double& theta, double& phi, double& chi)
{
  // get molecule coordinates WRT molecule center
  Array<int,1> mol1Members, mol2Members;
  PathData.Mol.MembersOf(mol1Members, mol1);
  PathData.Mol.MembersOf(mol2Members, mol2);
  dVec p1a = PathData.Path(slice1, mol1Members(1)) - PathData.Path(slice1, mol1Members(0));
  dVec p1b = PathData.Path(slice1, mol1Members(2)) - PathData.Path(slice1, mol1Members(0));
  dVec p2a = PathData.Path(slice2, mol2Members(1)) - PathData.Path(slice2, mol2Members(0));
  dVec p2b = PathData.Path(slice2, mol2Members(2)) - PathData.Path(slice2, mol2Members(0));
  // must impose minimum image for molecules at boundaries
  PathData.Path.PutInBox(p1a);
  PathData.Path.PutInBox(p1b);
  PathData.Path.PutInBox(p2a);
  PathData.Path.PutInBox(p2b);
	p1a = Normalize(p1a);
	p1b = Normalize(p1b);
	p2a = Normalize(p2a);
	p2b = Normalize(p2b);

  // compute bisector, norm, and in-plane unit vectors
  dVec b1 = GetBisector(p1a, p1b);
  dVec n1 = Normalize(cross(p1a, p1b));
  dVec u1 = Normalize(cross(b1, n1));
  dVec b2 = GetBisector(p2a, p2b);
  dVec n2 = Normalize(cross(p2a, p2b));
  dVec u2 = Normalize(cross(b2, n2));

  // theta is angle between in-plane vectors u
  // (rotation about norm to molecule plane)
  dVec A;
  if(isEqual(u1, u2)) {
    theta = 0.;
    phi = GetAngle(b2, b1);
    dVec N = Normalize(cross(b1, b2));
    if(isEqualOpp(N, u1))
      phi = 2*M_PI - phi;
    chi = 0.;
  }
  else if(isEqualOpp(u1, u2)) {
    theta = M_PI;
    phi = 0.;
    dVec newB = -1 * b1;
    chi = GetAngle(b2, newB);
    dVec N = Normalize(cross(newB, b2));
    if(isEqualOpp(N, u2))
      chi = 2*M_PI - chi;
  } else {
    theta = GetAngle(u2, u1);
    // axis of rotation for theta is cross of in-plane vectors u
    A = Normalize(cross(u1,u2));

    // phi is angle between norm n1 and A
    if(isEqual(A, n1))
      phi = 0.;
    else if(isEqualOpp(A, n1))
      phi = M_PI;
    else {
      phi = GetAngle(n1, A);
      dVec N = Normalize(cross(n1, A));
      if(isEqualOpp(N,u1))
        phi = 2*M_PI - phi;
    }
  
    // chi is angle between norm A and n2
    if(isEqual(A, n2))
      chi = 0.;
    else if(isEqualOpp(A, n2))
      chi = M_PI;
    else {
      chi = GetAngle(A, n2);
      dVec N = Normalize(cross(A, n2));
      if(isEqualOpp(N,u2))
        chi = 2*M_PI - chi;
    }
  }

  //cerr << "ANGLES " << phi << " " << theta << " " << chi << endl;
  //cerr << "       " << phiInv << " " << thetaInv << " " << chiInv << endl;

  //if(isnan(chi)) {
  //  cerr << "Crap chi=NAN A " << A << " " << n1 << " " << n2 << " theta " << theta <<  " " << endl;
  //  cerr << "isEqual gives " << isEqual(A,n2) << " and negative " << isEqualOpp(A, n2) << " between " << A << " " << (-1*n2) << endl;
  //  cerr << "p1a " << p1a << " p1b " << p1b << " p2a " << p2a << " p2b " << p2b << endl;
  //}
  //if(isnan(phi)) {
  //  cerr << "Crap phi=NAN A " << A << " " << n1 << " " << n2 << " theta " << theta <<  " " << endl;
  //  cerr << "isEqual gives " << isEqual(A,n1) << " and negative " << isEqualOpp(A, n1) << " between " << A << " " << (-1*n1) << endl;
  //  cerr << "p1a " << p1a << " p1b " << p1b << " p2a " << p2a << " p2b " << p2b << endl;
  //}
  if(phi<0){
    cerr << "NEGATIVE phi " << phi << endl;
    exit(1);
  }
  if(chi<0) {
    cerr << "NEGATIVE chi " << chi << endl;
    exit(1);
  }
}

FixedAxisRotorClass::FixedAxisRotorClass(PathDataClass &pathData ) : 
  RotorActionBaseClass (pathData)
{
  // hard-wired values
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  C = 1.;//0.5*2.65146*0.0493089;
}

double FixedAxisRotorClass::SingleAction (int slice1, int slice2, 
		       const Array<int,1> &changedParticles, int level)
{
  //cerr << "FA action intertia " << Ia << " " << Ib << " " << Ic << endl;
  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int numChangedMol = molRoster.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      //if(slice==1) {
      //  AltFixedAxis(slice, mol);
      //}
      
      double r, rAlt;
      rAlt = AltFixedAxis(slice, mol, phi, theta, chi);
      r = rAlt;
      //cerr << "GETANGLES_ALT " << theta << " " << phi << " " << chi << " rho " << rAlt << endl;

      //getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      //r = CalcRho(phi, theta, chi);
      //cerr << "GETANGLES " << theta << " " << phi << " " << chi << " rho " << r << endl;
      // need to check for negative rho; use noisy correction
      //if(r<0) {
      //  cout << "Got NEGATIVE rho " << r;
      //  double a = -1 * PathData.Path.Random.Local() * r;
      //  cout << " generated noise " << a << endl;
      //  r = a;
      //}

      double logRho = log(r);
      //cout << "Computed FIXED AXIS rho " << r << " and log " << logRho << endl;
      //double logRho = spline(theta, phi, chi);
      TotalK -= logRho;
    }
  }
  //cout << "returning FIXED AXIS Action " << TotalK << endl << endl;
  return (TotalK);
}

double FixedAxisRotorClass::d_dBeta (int slice1, int slice2, int level)
{
  int M = PathData.Path.NumTimeSlices()-1;
  //cerr << "Calculating FIXED AXIS energy over slices " << slice1 << " " << slice2 << " with " << M << " total slices" << endl;
  double TotalK = 0.0;
  double Z = 0.0;
  double altFA = 0.;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      //double r = CalcRho(phi, theta, chi);
      double FAR = AltFixedAxis (slice, mol, phi, theta, chi);
      TotalK += CalcFAEnergy(phi, theta, chi);
      //cout << theta << " " << phi << " " << chi << " ";// << endl;
      //altFA += AltFixedAxis(slice, mol);
      ////Z += r;
      //cout << altFA << " " << TotalK << endl;
    }
  }
  return (TotalK);
}

void FixedAxisRotorClass::Read (IOSectionClass &in)
{
  // read in moments of inertia
  // these are the defaults for water
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  in.ReadVar("Ia",Ia);
  in.ReadVar("Ib",Ib);
  in.ReadVar("Ic",Ic);

  // arbitrary constant (?)
  //C = 0.5*2.65146*0.0493089;
  C = 1.0;
}

string FixedAxisRotorClass::GetName()
{
  return "FixedAxisRotor";
}

double FixedAxisRotorClass::CalcRho(double phi, double theta, double psi)
{
  double tau = PathData.Path.tau;
  double chi = cos(theta/2) * cos((psi + phi)/2);
  double eta = sin(theta/2) * cos((psi - phi)/2);
  double xi = sin(theta/2) * sin((psi - phi)/2);
  double zeta = cos(theta/2) * sin((psi + phi)/2);

  double gamma = 2 * acos(chi);
  double sinG2 = sin(gamma/2);
  double gammaOverSinG2 = gamma/sinG2;
  //if(phi==0. && psi==0. && theta==0.)
  if(gamma<TOL) {
    gammaOverSinG2 = 2.;
  }

  double nDotI = (eta*eta*Ib + xi*xi*Ic + zeta*zeta*Ia);///(eta*eta + xi*xi + zeta*zeta);

  double SFA;
  SFA = 0.5 * gammaOverSinG2 * gammaOverSinG2 * nDotI/tau;
  //cerr << "FA ndoti " << nDotI << " tau " << tau << " gammaOSin " << gammaOverSinG2 << " SFA " << SFA << " angles " << phi << " " << theta << " " << psi << endl;

  double D = 0.5 * sqrt((Ia * Ib * Ic)/(tau*tau*tau))*gammaOverSinG2;

  double Q = D * exp(-SFA);
  //cerr << "C " << C << " D " << D << " SFA " << SFA << endl;
  return Q;
}

// testing
double FixedAxisRotorClass::AltFixedAxis (int slice, int mol, double& phi, double& theta, double& psi)
{
  double RotK = 0.0;
  int skip = 1;
  int ptcl1 = PathData.Mol.MembersOf(mol)(1);
  int ptcl2 = PathData.Mol.MembersOf(mol)(2);
  dVec P1 = PathData.Path(slice,ptcl1);
  dVec P2 = PathData.Path(slice,ptcl2);
  dVec P1prime = PathData.Path(slice+skip,ptcl1);
  dVec P2prime = PathData.Path(slice+skip,ptcl2);
  //cerr << "P1 " << P1 << endl;
  //cerr << "P2 " << P2 << endl;
  //cerr << "P1prime " << P1prime << endl;
  //cerr << "P2prime " << P2prime << endl;
  int Optcl = PathData.Mol(ptcl1);
  assert(mol == Optcl);
  dVec O = PathData.Path(slice,Optcl);
  dVec Oprime = PathData.Path(slice+skip,Optcl);
  //cerr << "identified COM " << Optcl << " with coords " << O << " and " << Oprime << endl;
  // Redefine coordinates WRT COM
  P1 -= O;
  P2 -= O;
  P1prime -= Oprime;
  P2prime -= Oprime;
  PathData.Path.PutInBox(P1);
  PathData.Path.PutInBox(P1prime);
  PathData.Path.PutInBox(P2);
  PathData.Path.PutInBox(P2prime);
  P1 = Normalize(P1);
  P2 = Normalize(P2);
  P1prime = Normalize(P1prime);
  P2prime = Normalize(P2prime);
  //cerr << "now WRT COM coords are " << endl;
  //cerr << "P1 " << P1 << endl;
  //cerr << "P2 " << P2 << endl;
  //cerr << "P1prime " << P1prime << endl;
  //cerr << "P2prime " << P2prime << endl;
  // Calculate bisectors for each configuration
  dVec n = Normalize(GetBisector(P1,P2));
  dVec nprime = Normalize(GetBisector(P1prime,P2prime));
  //cerr << "n " << n << endl;
  //cerr << "nprime " << nprime << endl;
  double vel_squared;
  double prefactor;
  //double theta, phi, psi;
  dVec r;
  if (isEqual(n, nprime)){
    cerr << "EQUAL--------------------------------" << endl;
    theta = 0.;
    r = Normalize(cross(P1prime,P2prime));
  }
  else{
// << "I DECIDED THEY WEREN'T EQUAL.  WHATEVER." << endl;
    // Calculate the cross product - the axis of rotation
    r = Normalize(cross(n,nprime));
// << "r " << r << endl;
    // Calculate polar angles and trig functions
    // Calculate azimuthal angle
    theta = GetAngle(n,nprime);
  }
  //cerr << "theta " << theta << endl;
  // Calculate lever arms and kinetic energy contributions (mass contained in lambda factor)
  //double alpha = HOH_half_angle;
  double alpha = 0.912135247; // 52.26deg
  double SinAlpha = sin(alpha);
  double CosAlpha = cos(alpha);
  dVec u1 = Normalize(P1 - Renormalize(n,CosAlpha));
  dVec z1 = Normalize(cross(P1,n));
  dVec u2 = Normalize(P2 - Renormalize(n,CosAlpha));
  dVec z2 = Normalize(cross(P2,n));
  dVec u1prime = Normalize(P1prime - Renormalize(nprime,CosAlpha));
  dVec z1prime = Normalize(cross(P1prime,nprime));
  dVec u2prime = Normalize(P2prime - Renormalize(nprime,CosAlpha));
  dVec z2prime = Normalize(cross(P2prime,nprime));
  double phi1 = GetAngle(z1,r);
  double phi2 = GetAngle(z2,r);
  double psi1 = GetAngle(z1prime,r);
  double psi2 = GetAngle(z2prime,r);
  //cerr << "P1 " << P1 << endl;
  //cerr << "n " << n << endl;
  //cerr << "P2 " << P2 << endl;
  phi = phi1;
  psi = psi1;
  //cerr << "ALT " << theta << " " << phi << " " << psi << " ";//<< endl;

  // Calculate quaternions
  double chi = cos(theta/2)*cos((psi + phi)/2);
  double eta = sin(theta/2)*cos((psi - phi)/2);
  double xi = sin(theta/2)*sin((psi - phi)/2);
  double zeta = cos(theta/2)*sin((psi + phi)/2);

  // Calculate Fixed-Axis Approximation parameters
  double gamma = 2*acos(chi);
  double SinHalfGamma = sin(gamma/2);
  // Rotation Moment Matrix
  //double Lx = R*R;
  //double Ly = R*R*CosAlpha*CosAlpha;
  //double Lz = R*R*SinAlpha*SinAlpha;
  double Lx = 0.027215;
  double Ly = 0.078226;
  double Lz = 0.04167;
  // Axis of Rotation;
  double ex = eta/SinHalfGamma;
  double ey = -xi/SinHalfGamma;
  double ez = zeta/SinHalfGamma;
  // Rotation Matrix Elements <e|L|e>
  double Mx = ex*Lx*ex;
  double My = ey*Ly*ey;
  double Mz = ez*Lz*ez;
  double MSum = Mx + My + Mz;
  vel_squared = MSum*gamma*gamma;
  prefactor = gamma/(2*SinHalfGamma);
  double GaussSum=0.0;
  RotK += prefactor*exp(-vel_squared);
  return (RotK);
}

double FixedAxisRotorClass::CalcFAEnergy(double phi, double theta, double psi)
{
  double tau = PathData.Path.tau;
  double chi = cos(theta/2) * cos((psi + phi)/2);
  double eta = sin(theta/2) * cos((psi - phi)/2);
  double xi = sin(theta/2) * sin((psi - phi)/2);
  double zeta = cos(theta/2) * sin((psi + phi)/2);

  double gamma = 2 * acos(chi);
  double sinG2 = sin(gamma/2);
  double gammaOverSinG2 = gamma/sinG2;
  //if(gammaOverSinG2 > 1e2) {
  //  cerr << phi << " " << theta << " " << psi << " " << gamma << " " << sinG2 << " " << gammaOverSinG2 << endl;
  //}
  //if((phi==0. && psi==0.) || theta==0.)
  //  gammaOverSinG2 = 2.;
  if(gamma<TOL) {
      gammaOverSinG2 = 2.;
  }

  double nDotI = (eta*eta*Ib + xi*xi*Ic + zeta*zeta*Ia);///(eta*eta + xi*xi + zeta*zeta);

  //double SFA;
  //SFA = 0.5 * gammaOverSinG2 * gammaOverSinG2 * nDotI/tau;

  //double D = 0.5 * sqrt((Ia * Ib * Ic)/(tau*tau*tau))*gammaOverSinG2;

  double Q = (-1.0 * nDotI * gammaOverSinG2 * gammaOverSinG2)/(2 * tau * tau)
    + 1.5/tau;

  return Q;
}
