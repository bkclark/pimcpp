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
  Array<double,1> &Vlong_r;
  double Z1Z2;
  int Level;
  inline double Vintegrand(double r)
  {
    double Vshort = Z1Z2/r - Vlong_r(r);
    return r*r*Vshort;
  }

public:
  inline double operator()(double r)
  {
    return Vintegrand(r);
  }
  UshortIntegrand (Array<double,1> &temp_Vlong_r,double &t_Z1Z2) :
    Vlong_r(temp_Vlong_r), Z1Z2(t_Z1Z2)
  { /* do nothing else */ }
};


class UlongIntegrand
{
private:
  Array<double,1> &Vlong_r;
  inline double Vintegrand(double r)
  {
    double Vlong = Vlong_r(r);
    return r*r*Vlong;
  }

public:
  inline double operator()(double r)
  {
      return Vintegrand(r);
  }
  UlongIntegrand (Array<double,1> &temp_Vlong_r) :
    Vlong_r(temp_Vlong_r)
  { /* do nothing else */ }
};


bool same(pair<double,double> &a, pair<double,double> &b)
{
  return fabs(a.first-b.first)<1e-6;
}

class EwaldClass
{
public:
  double Z1Z2;
  double kCut;
  dVec box;
  double vol;
  double r0;
  double rc;
  int NumImages;
  int NumPoints;
  int NumKnots;

  void SetkCut(double t_kCut)
  {
    kCut=t_kCut;
  }

  void SetCharge(double t_Z1Z2)
  {
    Z1Z2=t_Z1Z2;
  }

  void SetBox(dVec t_Box)
  {
    box=t_Box;
    vol=1.0;
    for (int i=0; i<DIM; i++)
      vol *= box[i];
  }

  void Setr0(double t_r0)
  {
    r0=t_r0;
    if (r0 == 0.0)
      cerr<<"Warning: r0 = 0.0!"<<endl;
  }

  void SetrCut(double t_rc)
  {
    rc=t_rc;
  }

  void SetNumPoints(double t_NumPoints)
  {
    NumPoints=t_NumPoints;
  }

  void SetNumImages(double t_NumImages)
  {
    NumImages=t_NumImages;
  }

  void SetNumKnots(double t_NumKnots)
  {
    NumKnots=t_NumKnots;
  }

  double Xk_V (double k,double rcut)
  {
    double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
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

  void SetupkVecs(Array<dVec,1> &kVecs, Array<double,1> &MagK, bool includeAll=false)
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
    cerr << "kCutoff = " << kCutoff << endl;
    cerr << "Number of kVecs = " << numVecs << endl;
    cerr << "MaxkIndex = " << MaxkIndex << endl;
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

  double Breakup()
  {
    const double tolerance = 1.0e-11;
    double boxVol = box[0]*box[1]*box[2];
    //for (int i=1; i<DIM; i++)
    //  rc = min (rc, 0.5*box[i]);
    double kvol = 1.0; // Path.GetkBox()[0];
    for (int i=0; i<DIM; i++)
      kvol *= 2*M_PI/box[i]; // Path.GetkBox()[i];
    double kavg = pow(kvol,1.0/3.0);

    LPQHI_BasisClass basis;
    basis.Set_rc(rc);
    basis.SetBox(box);
    basis.SetNumKnots(NumKnots);

    // We try to pick kcont to keep reasonable number of k-vectors
    double kCont = 50.0 * kavg;
    double delta = basis.GetDelta();
    double kMax = 20.0*M_PI/delta;

    OptimizedBreakupClass breakup(basis);
    breakup.SetkVecs (kCut, kCont, kMax);
    int numk = breakup.kpoints.size();
    int N = basis.NumElements();
    Array<double,1> t(N);
    Array<bool,1>   adjust (N);
    Array<double,1> Xk(numk);
    Array<double,1> Vlong_k;

    // Would be 0.5, but with two timeslice distdisp, it could be a
    // little longer
    double rmax = 0.75 * sqrt (dot(box,box));
    LinearGrid LongGrid;
    LongGrid.Init (r0, rmax, NumPoints);
    Array<double,1> Vlong_r(NumPoints), r(NumPoints), Vshort_r(NumPoints);
    for (int i=0; i<NumPoints; i++)
      r(i) = LongGrid(i);
    Array<dVec,1> kVecs;
    Array<double,1> kMag;
    SetupkVecs(kVecs,kMag,false);
    Vlong_k.resize(kVecs.size()); //need actual number of ks
    Vlong_k=0.0;
    Vlong_r = 0.0;

    // Calculate Xk's
    for (int ki=0; ki<numk; ki++) {
      // Xk(ki) = CalcXk(paIndex, 0, breakup.kpoints(ki)[0], rc, JOB_V);
      double k = breakup.kpoints(ki)[0];
      Xk(ki) = Xk_V (k,rc) / boxVol;
    }

    // Set boundary conditions at rc:  Force value and first and
    // second derivatives of long-range potential to match the full
    // potential at rc.
    adjust = true;
    delta = basis.GetDelta();
    //     t(N-3) = pa.V(rc);                 adjust(N-3) = false;
    //     t(N-2) = pa.Vp(rc)*delta;          adjust(N-2) = false;
    //     t(N-1) = pa.Vpp(rc)*delta*delta;   adjust(N-1) = false;
    //     t(1) = 0.0;                        adjust(1)   = false;

    //     t(N-3) = 1.0/delta;                 adjust(N-3) = false;
    //     t(N-2) = (-1.0/(delta*delta))*delta;          adjust(N-2) = false;
    //     t(N-1) = (2.0/(delta*delta*delta))*delta*delta;   adjust(N-1) = false;
    //     t(1) = 0.0;                        adjust(1)   = false;


    // Now, do the optimal breakup:  this gives me the coefficents
    // of the basis functions, h_n in the array t.
    cerr << "Doing V breakup...\n";
    double chi=breakup.DoBreakup (Xk, t, adjust);
    cerr<<"Chi is "<<chi<<endl;
    cerr << "Done.\n";

    // Now, we must put this information into the pair action
    // object.  First do real space part
    double Vlong_r0=0.0;
    for (int n=0; n<N; n++)
      Vlong_r0 += t(n)*basis.h(n,0.0);
    for (int i=0; i<LongGrid.NumPoints; i++) {
       double r = LongGrid(i);
       if (r <= rc) {
         // Sum over basis functions
         for (int n=0; n<N; n++)
           Vlong_r(i) += t(n) * basis.h(n, r);
       }
       else
         Vlong_r(i) = Z1Z2/r; // 0.0; //hack!  pa.V (r);
    }

    // Calculate FT of Ushort at k=0
    UshortIntegrand shortIntegrand(Vlong_r,Z1Z2);
    GKIntegration<UshortIntegrand, GK31> shortIntegrator(shortIntegrand);
    shortIntegrator.SetRelativeErrorMode();
    double Vshort_k0 = 4.0*M_PI/boxVol * shortIntegrator.Integrate(1.0e-100, rc, tolerance);
    cerr << "Vshort_k0 = " << Vshort_k0 << endl;

    // Calculate FT of Vlong at k=0
    UlongIntegrand longIntegrand(Vlong_r);
    GKIntegration<UlongIntegrand, GK31> longIntegrator(longIntegrand);
    longIntegrator.SetRelativeErrorMode();
    double Vlong_k0 = 4.0*M_PI/boxVol * longIntegrator.Integrate(1.0e-100, rmax, tolerance);
    cerr << "Vlong_k0 = " << Vlong_k0 << endl;

    ofstream outfile;
    outfile.open("rData.txt");
    for (int i=0;i<LongGrid.NumPoints;i++){
      double r=LongGrid(i);
      outfile<<r<<" "<<Z1Z2/r-Vlong_r(i)<<" "<<Vlong_r(i)<<endl;
    }
    outfile.close();

    // Now do k-space part
    outfile.open("kData.txt");
    outfile<<0.0<<" "<<(Vlong_k0+Vshort_k0)<<endl;
    for (int ki=0; ki < kVecs.size(); ki++) {
      double k = sqrt(dot(kVecs(ki),kVecs(ki)));
      // Sum over basis functions
      for (int n=0; n<N; n++)
        Vlong_k(ki) += t(n) * basis.c(n,k);
      // Now add on part from rc to infinity
      //pa.Vlong_k(ki) -= CalcXk(paIndex, 0, k, rc, JOB_V);
      Vlong_k(ki) -= Xk_V(k,rc) / boxVol;
      outfile<<k<<" "<<Vlong_k(ki)*boxVol<<endl;
    }
    outfile.close();
  }

  void TestMadelung()
  {
    // Get Long-ranged k-space part
    cout << "Getting long-ranged k-space part" << endl;
    ifstream infile;
    infile.open("kData.txt");
    vector<double> ks;
    vector<double> fVls;
    double Vk0=0.0;
    while (!infile.eof()) {
      double k;
      double fVl;
      infile>>k;
      if (!infile.eof()) {
        infile>>fVl;
        ks.push_back(k);
        fVls.push_back(fVl/Z1Z2);
        if (k == 0.0)
          Vk0 = fVl/Z1Z2;
      }
    }

    // Get Short-ranged r-space part
    cout << "Getting short-ranged r-space part" << endl;
    ifstream pointFile;
    pointFile.open("rData.txt");
    Array<double,1> rs(NumPoints);
    Array<double,1> Vss(NumPoints);
    int ri = 0;
    while (!pointFile.eof()) {
      double r;
      double Vs;
      pointFile>>r;
      pointFile>>Vs;
      rs(ri) = r;
      Vss(ri) = Vs/Z1Z2;
      ri += 1;
    }
    double r0 = rs(0);
    double rMax = rs(rs.size()-1);
    cout << "Creating spline" << endl;
    LinearGrid LongGrid;
    LongGrid.Init(r0, rMax, NumPoints);
    CubicSpline VsSpline(&LongGrid,Vss);

    cout << "Setting coordinates" << endl;
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
    cout << "Getting Vs" << endl;
    double L = box[0];
    double Li = 1.0/L;
    double Vs = 0.;
    for (int i=0; i<xs.size()-1; i++) {
      for (int j=i+1; j<xs.size(); j++) {
        dVec r = xs(i) - xs(j);
        for (int d=0; d<DIM; d++)
          r(d) -= floor(r(d)*Li + 0.5)*L; // Put in box
        double magr = sqrt(dot(r,r));
        if (magr <= rc)
          Vs += Qs(i)*Qs(j)*VsSpline(magr);
      }
    }
    cout << "Vs " << Vs << endl;

    // Long-ranged k-space
    Array<dVec,1> kVecs;
    Array<double,1> MagK;
    SetupkVecs(kVecs,MagK,true);
    Array<double,1> fVl(kVecs.size());
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
    cout << "Vl " << Vl << endl;

    // Self-interacting terms
    double Vself = 0;
    double alpha = 7./box[0]; // HACK !! ! ! ! ! !!  ! !! ! ! ! ! ! ! ! !!!!
    for (int i=0; i<Qs.size(); i++)
      Vself -= Qs(i)*Qs(i)*alpha/sqrt(M_PI);
    cout << "Vself " << Vself << endl;

    // Neutralizing Background
    double rMin = 1e-5;
    rMax = rc;
    int nR = 10000;
    double dr = (rMax-rMin)/nR;
    LinearGrid nbGrid;
    nbGrid.Init(rMin, rc, nR);
    double s0 = 0.;
    for (int i=0; i<nbGrid.NumPoints; i++) {
      double r = nbGrid(i);
      s0 += dr*4*M_PI*r*r*VsSpline(r)/vol;
    }

    double Vb = 0.;
    for (int i=0; i<N; i++)
      Vb -= 0.5*N*N*s0;
    cout << "Vb " << Vb << endl;

    double V = Vs + Vl + Vself;

    // Madelung Constant
    double estMad = box[0]*V/N;
    cout << "Estimated madelung constant: " << estMad << endl;
    double extMad = -1.747564594633182190636212035544397403481;
    cout << "Exact madelung constant: " << extMad << endl;
    cout << "Relative error: " << abs(estMad-extMad)/extMad << endl;
  }

  void BasicEwald(double &alpha)
  {
    Array<dVec,1> kVecs;
    Array<double,1> MagK;
    SetupkVecs(kVecs,MagK,true);
    ofstream outfile;

    // Short-ranged r-space part
    outfile.open("rData.txt");
    double rMax = rc;
    Array<double,1> Vss(NumPoints);
    LinearGrid LongGrid;
    LongGrid.Init (r0, rMax, NumPoints);
    for (int i=0; i<LongGrid.NumPoints; i++){
      double r = LongGrid(i);
      Vss(i) = Z1Z2*erfc(alpha*r)/r;
      outfile<<r<<" "<<Vss(i)<<endl;
    }
    outfile.close();

    // Long-ranged k-space part
    outfile.open("kData.txt");
    double Vlong_k0 = -4.0*M_PI*Z1Z2/(4.0*alpha*alpha*vol);
    outfile<<0.0<<" "<<Vlong_k0<<endl;
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

};

int main(int argc, char* argv[])
{
  EwaldClass e;
  double L = atof(argv[1]);
  dVec box(L,L,L);
  e.SetBox(box);
  e.SetkCut(2.0*M_PI*atof(argv[2])/L);
  e.Setr0(atof(argv[3]));
  e.SetrCut(atof(argv[4]));
  e.SetCharge(atof(argv[5]));
  e.SetNumPoints(atoi(argv[6]));
  bool doBreakup = atoi(argv[7]);

  if (doBreakup) {
    e.SetNumKnots(atoi(argv[8]));
    e.Breakup();
    e.SetNumImages(0);
  } else {
    double alpha=7.0/(L);
    e.SetNumImages(50);
    e.BasicEwald(alpha);
  }

  e.TestMadelung();
}
