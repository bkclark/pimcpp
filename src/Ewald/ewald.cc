#include "../Blitz.h"
#include "OptimizedBreakup.h"
#include "../Splines/Grid.h"
#include "../Splines/CubicSpline.h"
#include <vector>
#include <algorithm>
#include "../Integration/GKIntegration.h"
#include <complex>

const int DIM = 3;

class UshortIntegrand
{
private:
  Array<double,1> &Vl;
  double Z1Z2;
  int Level;
  inline double Vintegrand(double r)
  {
    return r*r*(Z1Z2/r - Vl(r));
  }

public:
  inline double operator()(double r)
  {
    return Vintegrand(r);
  }
  UshortIntegrand (Array<double,1> &temp_Vl,double &t_Z1Z2) :
    Vl(temp_Vl), Z1Z2(t_Z1Z2)
  { /* do nothing else */ }
};


class UlongIntegrand
{
private:
  Array<double,1> &Vl;
  inline double Vintegrand(double r)
  {
    return r*r*Vl(r);
  }

public:
  inline double operator()(double r)
  {
      return Vintegrand(r);
  }
  UlongIntegrand (Array<double,1> &temp_Vl) :
    Vl(temp_Vl)
  { /* do nothing else */ }
};


bool same(pair<double,double> &a, pair<double,double> &b)
{
  return fabs(a.first-b.first)<1e-6;
}

class EwaldClass
{
public:
  dVec box;
  double Z1Z2, vol, rMin, rMax;
  Array<dVec,1> kVecs;
  Array<double,1> MagK;
  double kCut;
  int nMax, nPoints, nKnots;
  int gridType, breakupType;
  Grid &grid;

  EwaldClass(double t_Z1Z2, dVec t_Box, int t_nMax, double t_rMin, double t_rMax, int t_nPoints,
             Grid &t_grid, int t_breakupType, int t_nKnots=0)
    : Z1Z2(t_Z1Z2), box(t_Box), nMax(t_nMax), rMin(t_rMin), rMax(t_rMax), nPoints(t_nPoints),
      grid(t_grid), breakupType(t_breakupType), nKnots(t_nKnots)
  {
    if (rMin == 0.)
      cerr << "Warning: rMin = 0!" << endl;
    vol = 1.;
    for (int i=0; i<DIM; i++)
      vol *= box[i];
    kCut = 2.*M_PI*nMax/box[0];
    SetupkVecs();
  }

  double DoBreakup()
  {
    if (breakupType == 0)
      BasicEwald(7./box[0]); // HACK: Hard-coded alpha
    else if (breakupType == 1)
      OptimizedBreakup();
    else
      cerr << "ERROR: Unrecognized breakup type." << endl;
  }

  double Xk_V (double k,double rMaxut)
  {
    double C0 = -4.0*M_PI/(k*k) * cos(k*rMaxut);
    return (Z1Z2*C0);
  }

  inline bool Include(dVec k, bool includeAll=false)
  {
    if (includeAll)
      return true;
    //  assert (DIM == 3);
    if (abs(k[0]) > 0.0)
      return true;
    else if ((k[0]==0.0) && (abs(k[1]) > 0.0))
      return true;
    else if ((DIM==3) && ((k[0]==0.0) && (k[1]==0.0) && (abs(k[2]) > 0.0)))
      return true;
    else
      return false;
  }

  void SetupkVecs(bool includeAll=false)
  {
    double kCutoff=kCut;
    //    Array<double,1> MagK;
    dVec kBox;
    for (int i=0; i<DIM; i++)
      kBox[i] = 2.0*M_PI/box[i];

    int numVecs=0;
    dVec MaxkIndex;
    for (int i=0;i<DIM;i++){
      MaxkIndex[i]= (int) ceil(1.*kCutoff/kBox[i]);
    }

    dVec k;
    TinyVector<int,DIM> ki;
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      k[0] = ix*kBox[0];
      if (DIM > 1) {
        for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
          k[1] = iy*kBox[1];
          if (DIM > 2) {
            for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
              k[2] = iz*kBox[2];
              if (Include(k,includeAll))
                numVecs++;
            }
          } else {
            if (Include(k,includeAll))
              numVecs++;
          }
        }
      } else {
        if (Include(k,includeAll))
           numVecs++;
      }
    }
    //  kIndices.resize(numVecs);
    cerr << "kCutoff = " << kCutoff << ", # of kVecs = " << numVecs << endl;
    kVecs.resize(numVecs);
    MagK.resize(numVecs);

    //   for (int i=0; i<DIM; i++)
    //     C[i].resize(2*MaxkIndex[i]+1);
    numVecs = 0;
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      k[0] = ix*kBox[0];
      ki[0]= ix+MaxkIndex[0];
      if (DIM > 1) {
        for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
          k[1] = iy*kBox[1];
          ki[1]= iy+MaxkIndex[1];
          if (DIM > 2) {
            for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
              k[2] = iz*kBox[2];
              ki[2]= iz+MaxkIndex[2];
              if (Include(k,includeAll)) {
                kVecs(numVecs) = k;
                MagK(numVecs) = sqrt(dot(k,k));
                numVecs++;
              }
            }
          } else {
            if (Include(k,includeAll)) {
              kVecs(numVecs) = k;
              MagK(numVecs) = sqrt(dot(k,k));
              numVecs++;
            }
          }
        }
      } else {
        if (Include(k,includeAll)) {
          kVecs(numVecs) = k;
          MagK(numVecs) = sqrt(dot(k,k));
          numVecs++;
        }
      }
    }
  }

  double OptimizedBreakup()
  {
    const double tolerance = 1.0e-11;
    double boxVol = box[0]*box[1]*box[2];
    double kvol = 1.0; // Path.GetkBox()[0];
    for (int i=0; i<DIM; i++)
      kvol *= 2*M_PI/box[i]; // Path.GetkBox()[i];
    double kavg = pow(kvol,1.0/3.0);

    LPQHI_BasisClass basis;
    basis.Set_rc(rMax);
    basis.SetBox(box);
    basis.SetNumKnots(nKnots);

    // We try to pick kcont to keep reasonable number of k-vectors
    double kCont = 50.0 * kavg;
    double delta = basis.GetDelta();
    double kMax = 20.0*M_PI/delta;

    OptimizedBreakupClass breakup(basis);
    breakup.SetkVecs (kCut, kCont, kMax);
    int numk = breakup.kpoints.size();
    int N = basis.NumElements();
    Array<double,1> t(N);
    Array<bool,1> adjust (N);
    Array<double,1> Xk(numk);
    Array<double,1> fVl;

    // Would be 0.5, but with two timeslice distdisp, it could be a
    // little longer
    Array<double,1> Vl(nPoints), r(nPoints), Vs(nPoints);
    for (int i=0; i<nPoints; i++)
      r(i) = grid(i);
    fVl.resize(kVecs.size()); //need actual number of ks
    fVl = 0.0;
    Vl = 0.0;

    // Calculate Xk's
    for (int ki=0; ki<numk; ki++) {
      // Xk(ki) = CalcXk(paIndex, 0, breakup.kpoints(ki)[0], rMax, JOB_V);
      double k = breakup.kpoints(ki)[0];
      Xk(ki) = Xk_V (k,rMax) / boxVol;
    }

    // Set boundary conditions at rMax:  ForMaxe value and first and
    // second derivatives of long-range potential to match the full
    // potential at rMax.
    adjust = true;
    delta = basis.GetDelta();
    //     t(N-3) = pa.V(rMax);                 adjust(N-3) = false;
    //     t(N-2) = pa.Vp(rMax)*delta;          adjust(N-2) = false;
    //     t(N-1) = pa.Vpp(rMax)*delta*delta;   adjust(N-1) = false;
    //     t(1) = 0.0;                        adjust(1)   = false;

    //     t(N-3) = 1.0/delta;                 adjust(N-3) = false;
    //     t(N-2) = (-1.0/(delta*delta))*delta;          adjust(N-2) = false;
    //     t(N-1) = (2.0/(delta*delta*delta))*delta*delta;   adjust(N-1) = false;
    //     t(1) = 0.0;                        adjust(1)   = false;


    // Now, do the optimal breakup:  this gives me the coefficents
    // of the basis functions, h_n in the array t.
    double chi = breakup.DoBreakup (Xk, t, adjust);
    cerr<<"Chi = "<<chi<<endl;

    // Now, we must put this information into the pair action
    // object.  First do real space part
    double Vl0=0.0;
    for (int n=0; n<N; n++)
      Vl0 += t(n)*basis.h(n,0.0);
    cout << "Vl0 = " << Vl0 << endl;
    for (int i=0; i<grid.NumPoints; i++) {
       double r = grid(i);
       if (r <= rMax) {
         // Sum over basis functions
         for (int n=0; n<N; n++)
           Vl(i) += t(n) * basis.h(n, r);
       }
       else
         Vl(i) = Z1Z2/r; // 0.0; //hack!  pa.V (r);
    }

    // Calculate FT of Ushort at k=0
    UshortIntegrand shortIntegrand(Vl,Z1Z2);
    GKIntegration<UshortIntegrand, GK31> shortIntegrator(shortIntegrand);
    shortIntegrator.SetRelativeErrorMode();
    double fVs0 = 4.0*M_PI/boxVol * shortIntegrator.Integrate(1.0e-100, rMax, tolerance);
    cerr << "fVs0 = " << fVs0 << endl;

    // Calculate FT of Vl at k=0
    UlongIntegrand longIntegrand(Vl);
    GKIntegration<UlongIntegrand, GK31> longIntegrator(longIntegrand);
    longIntegrator.SetRelativeErrorMode();
    double fVl0 = 4.0*M_PI/boxVol * longIntegrator.Integrate(1.0e-100, rMax, tolerance);
    cerr << "fVl0 = " << fVl0 << endl;

    ofstream outfile;
    outfile.open("rData.txt");
    outfile<<0.<<" "<<Vl0<<endl;
    for (int i=0;i<grid.NumPoints; i++){
      double r=grid(i);
      outfile<<r<<" "<<Z1Z2/r-Vl(i)<<endl;
    }
    outfile.close();

    // Now do k-space part
    outfile.open("kData.txt");
    outfile<<0.0<<" "<<(fVl0+fVs0)<<endl;
    vector<pair<double, double> > fVls;
    for (int ki=0; ki < kVecs.size(); ki++) {
      double k = sqrt(dot(kVecs(ki),kVecs(ki)));
      // Sum over basis functions
      for (int n=0; n<N; n++)
        fVl(ki) += t(n) * basis.c(n,k);
      // Now add on part from rMax to infinity
      //pa.fVl(ki) -= CalcXk(paIndex, 0, k, rMax, JOB_V);
      fVl(ki) -= Xk_V(k,rMax) / boxVol;
      fVls.push_back(make_pair(k, fVl(ki)));
    }
    sort(fVls.begin(), fVls.end());
    vector<pair<double, double> >::iterator new_end = unique(fVls.begin(), fVls.end(), same);
    // delete all elements past new_end
    fVls.erase(new_end, fVls.end());
    for (int i=1;i<fVls.size();i++){
      outfile<<fVls[i].first<<" "<<fVls[i].second<<endl; // *vol<<endl;
    }
    outfile.close();
  }


  void BasicEwald(double alpha)
  {

    // Short-ranged r-space part
    ofstream outfile;
    outfile.open("rData.txt");
    double Vl0 = Z1Z2*2.*alpha/sqrt(M_PI);
    outfile<<0.<<" "<<Vl0<<endl;
    Array<double,1> Vss(nPoints);
    for (int i=0; i<grid.NumPoints; i++){
      double r = grid(i);
      Vss(i) = Z1Z2*erfc(alpha*r)/r;
      outfile<<r<<" "<<Vss(i)<<endl;
    }
    outfile.close();

    // Long-ranged k-space part
    outfile.open("kData.txt");
    double fVl0 = -4.0*M_PI*Z1Z2/(4.0*alpha*alpha*vol);
    outfile<<0.0<<" "<<fVl0<<endl;
    vector<pair<double, double> > fVls;
    for (int i=0; i<kVecs.size(); ++i) {
      dVec k = kVecs(i);
      double k2 = dot(k,k);
      double kMag = sqrt(k2);
      double val = (4.*M_PI*Z1Z2/(k2*vol)) * exp(-k2/(4.*alpha*alpha));
      fVls.push_back(make_pair(kMag, val));
    }
    sort(fVls.begin(), fVls.end());
    vector<pair<double, double> >::iterator new_end = unique(fVls.begin(), fVls.end(), same);
    // delete all elements past new_end
    fVls.erase(new_end, fVls.end());
    for (int i=1;i<fVls.size();i++){
      outfile<<fVls[i].first<<" "<<fVls[i].second<<endl; // *vol<<endl;
    }
    outfile.close();

  }


  void TestMadelung()
  {
    // Get Long-ranged k-space part
    ifstream infile;
    infile.open("kData.txt");
    vector<double> ks;
    vector<double> fVls;
    double fV0 = 0.;
    while (!infile.eof()) {
      double k;
      double fVl;
      infile>>k;
      if (!infile.eof()) {
        infile>>fVl;
        ks.push_back(k);
        fVls.push_back(fVl/Z1Z2);
        if (k == 0.0)
          fV0 = fVl/Z1Z2;
      }
    }
    infile.close();

    // Get Short-ranged r-space part
    ifstream pointFile;
    pointFile.open("rData.txt");
    Array<double,1> rs(nPoints);
    Array<double,1> Vss(nPoints);
    int ri = 0;
    double Vl0 = 0.;
    while (!pointFile.eof()) {
      double r;
      double Vs;
      pointFile>>r;
      if (!pointFile.eof()) {
        pointFile>>Vs;
        if (r == 0.) {
          Vl0 = Vs/Z1Z2;
        } else {
          rs(ri) = r;
          Vss(ri) = Vs/Z1Z2;
          ri += 1;
        }
      }
    }
    pointFile.close();
    CubicSpline VsSpline(&grid,Vss);

    // Set Coordinates
    double s1 = box[0]/2.;
    int N = 8;
    Array<dVec,1> xs(N);
    Array<double,1> Qs(N);
    if (DIM==3) {
      dVec x0(0,0,0);
      xs(0) = x0;
      dVec x1(s1,s1,0);
      xs(1) = x1;
      dVec x2(s1, 0, s1);
      xs(2) = x2;
      dVec x3(0, s1, s1);
      xs(3) = x3;
      dVec x4(s1, 0, 0);
      xs(4) = x4;
      dVec x5(0, s1, 0);
      xs(5) = x5;
      dVec x6(0, 0, s1);
      xs(6) = x6;
      dVec x7(s1, s1, s1);
      xs(7) = x7;
      Qs(0) = 1.0;
      Qs(1) = 1.0;
      Qs(2) = 1.0;
      Qs(3) = 1.0;
      Qs(4) = -1.0;
      Qs(5) = -1.0;
      Qs(6) = -1.0;
      Qs(7) = -1.0;
    }

    // Short-ranged r-space
    double L = box[0];
    double Li = 1.0/L;
    double Vs = 0.;
    for (int i=0; i<xs.size()-1; i++) {
      for (int j=i+1; j<xs.size(); j++) {
        dVec r = xs(i) - xs(j);
        for (int d=0; d<DIM; d++)
          r(d) -= floor(r(d)*Li + 0.5)*L; // Put in box
        double magr = sqrt(dot(r,r));
        if (magr <= rMax)
          Vs += Qs(i)*Qs(j)*VsSpline(magr);
      }
    }
    cout << "Vs = " << Vs << endl;

    // Long-ranged k-space
    Array<double,1> fVl(kVecs.size());
    fVl = 0.;
    for (int i=0; i<kVecs.size(); i++) {
      bool foundMe = false;
      for (int j=0; j<ks.size(); j++) {
        if (fabs(MagK(i)-ks[j])<1.e-4) {
          if (!foundMe) {
            foundMe = true;
            fVl(i) = fVls[j];
          }
        }
      }
    }

    double Vl = 0.;
    for (int i=0; i<kVecs.size(); i++) {
      if(MagK(i) < kCut) {
        double Re = 0.;
        double Im = 0.;
        for (int j=0; j<xs.size(); j++) {
          double h = dot(kVecs(i), xs(j));
          Re += Qs(j)*cos(h);
          Im -= Qs(j)*sin(h);
        }
        for (int j=0; j<xs.size(); j++) {
          double h = dot(kVecs(i), xs(j));
          Vl += Qs(j) * (Re*cos(h) - Im*sin(h)) * fVl(i);
        }
      }
    }
    Vl *= 0.5;
    cout << "Vl = " << Vl << endl;

    // Self-interacting terms
    double Vself = 0;
    for (int i=0; i<Qs.size(); i++)
      Vself -= Qs(i)*Qs(i)*0.5*Vl0;
    cout << "Vself = " << Vself << endl;

    // Neutralizing Background
    double Vb = 0.;
    for (int i=0; i<N; i++)
      Vb += 0.5*N*N*fV0;
    cout << "Vbackground = " << Vb << endl;

    double V = Vs + Vl + Vself;
    cout << "V = Vs + Vl + Vself = " << V << endl;

    // Madelung Constant
    double estMad = box[0]*V/N;
    printf("Estimated madelung constant = %0.15f\n", estMad);
    double extMad = -1.747564594633182190636212035544397403481;
    printf("Exact madelung constant = %0.15f\n", extMad);
    cout << "Absolute error = " << abs(estMad-extMad) << endl;
    cout << "Relative error = " << abs(estMad-extMad)/extMad << endl;
  }

};

int main(int argc, char* argv[])
{
  double L = atof(argv[1]);
  dVec box(L,L,L); // Hard-coded cube
  double nMax = atoi(argv[2]);
  double rMin = atof(argv[3]);
  double rMax = atof(argv[4]);
  int nPoints = atoi(argv[5]);
  int gridType = atoi(argv[6]);

  Grid *grid;
  if (gridType == 0) {
    grid = new LinearGrid(rMin, rMax, nPoints);
  } else if (gridType == 1) {
    grid = new LogGrid(rMin, pow(rMax/rMin,1./nPoints), nPoints);
  }

  double Z1Z2 = atof(argv[7]);
  bool breakupType = atoi(argv[8]);
  int nKnots = atoi(argv[9]);

  EwaldClass e(Z1Z2, box, nMax, rMin, rMax, nPoints, *grid, breakupType, nKnots);
  e.DoBreakup();
  e.TestMadelung();

}
