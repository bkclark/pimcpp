#include "../Blitz.h"
#include "OptimizedBreakup.h"
#include "../Splines/Grid.h"
#include <vector>
#include <algorithm>
#include "../Integration/GKIntegration.h"


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
  double r0;
  double rc;
  int NumImages;
  int NumPoints;

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

  double Xk_V (double k,double rcut)
  {
    double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
    return (Z1Z2*C0);
  }

  inline bool Include(dVec k,bool includeAll=false)
  {
    if (includeAll)
      return true;
    //  assert (NDIM == 3);
    if (k[0] > 0.0)
      return true;
    else if ((k[0]==0.0) && (k[1]>0.0))
      return true;
    else if ((NDIM==3) && ((k[0]==0.0) && (k[1]==0.0) && (k[2] > 0.0)))
      return true;
    else
    return false;
  }

  void SetupkVecs3D(Array<dVec,1> &kVecs, Array<double,1> &MagK, bool includeAll=false)
  {
    double kCutoff=kCut;
    assert (NDIM == 3);
    //    Array<double,1> MagK;
    dVec kBox;
    for (int i=0; i<NDIM; i++)
      kBox[i] = 2.0*M_PI/box[i];

    int numVecs=0;
    dVec MaxkIndex;
    for (int i=0;i<NDIM;i++){
      MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    }

    dVec k;
    TinyVector<int,NDIM> ki;
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      k[0] = ix*kBox[0];
      for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
        k[1] = iy*kBox[1];
        for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
          k[2] = iz*kBox[2];
          if ((dot(k,k)<kCutoff*kCutoff) && Include(k,includeAll))
            numVecs++;
        }
      }
    }
    //  kIndices.resize(numVecs);
    cerr << "kCutoff = " << kCutoff << endl;
    cerr << "Number of kVecs = " << numVecs << endl;
    cerr << "MaxkIndex = " << MaxkIndex << endl;
    kVecs.resize(numVecs);
    MagK.resize(numVecs);

    //   for (int i=0; i<NDIM; i++)
    //     C[i].resize(2*MaxkIndex[i]+1);
    numVecs = 0;
    for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
      k[0] = ix*kBox[0];
      ki[0]= ix+MaxkIndex[0];
      for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
        k[1] = iy*kBox[1];
        ki[1]= iy+MaxkIndex[1];
        for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
          k[2] = iz*kBox[2];
          ki[2]= iz+MaxkIndex[2];
          if ((dot(k,k)<kCutoff*kCutoff) && Include(k,includeAll)) {
            kVecs(numVecs) = k;
            //          kIndices(numVecs)=ki;
            MagK(numVecs) = sqrt(dot(k,k));
            numVecs++;
          }
        }
      }
    }
  }

  double Breakup()
  {
    const double tolerance = 1.0e-11;
    double boxVol = box[0]*box[1]*box[2];
    for (int i=1; i<NDIM; i++)
      rc = min (rc, 0.5*box[i]);
    double kvol = 2*M_PI/box[0]; // Path.GetkBox()[0];
    for (int i=1; i<NDIM; i++)
      kvol *= 2*M_PI/box[i]; // Path.GetkBox()[i];
    double kavg = pow(kvol,1.0/3.0);
    int numKnots=8;

    LPQHI_BasisClass basis;
    basis.Set_rc(rc);
    basis.SetBox(box);
    basis.SetNumKnots (numKnots);

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
    SetupkVecs3D(kVecs,kMag,false);
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
    outfile<<0.0<<" "<<(Vlong_k0+Vshort_k0)/boxVol<<endl;
    for (int ki=0; ki < kVecs.size(); ki++) {
      const dVec &kv = kVecs(ki);
      double k = sqrt (dot(kv,kv));
      // Sum over basis functions
      for (int n=0; n<N; n++)
      Vlong_k(ki) += t(n) * basis.c(n,k);
      // Now add on part from rc to infinity
      //pa.Vlong_k(ki) -= CalcXk(paIndex, 0, k, rc, JOB_V);
      Vlong_k(ki) -= Xk_V(k,rc) / boxVol;
      outfile<<sqrt(kv[0]*kv[0]+kv[1]*kv[1]+kv[2]*kv[2])<<" ";
      outfile<<Vlong_k(ki)<<endl;
    }
    outfile.close();
  }

  void FourierTransform(bool diagonal=false)
  {
    ifstream infile;
    infile.open("kData.txt");
    vector<double> magDavid;
    vector<double> uk;
    double Vk0=0.0;
    while (!infile.eof()){
      double k;
      double data;
      infile>>k;
      if (!infile.eof()){
        infile>>data;
        magDavid.push_back(k);
        uk.push_back(data*(box[0]*box[1]*box[2]));
        if (k==0.0)
          Vk0 = data;
      }
    }
    Array<dVec,1> kVecs;
    Array<double,1> MagK;
    SetupkVecs3D(kVecs,MagK,true);
    ifstream pointFile;
    pointFile.open("rData.txt");
    ofstream outFile;
    outFile.open("FTData.txt");
    while (!pointFile.eof()){
      double r;
      double VShort_r;
      pointFile>>r;
      pointFile>>VShort_r;
      double VLong_r;
      pointFile>>VLong_r;
      double V=0.0;
      for (int i=0;i<kVecs.size();i++){

        bool foundMe=false;
        if ((MagK(i)==0))
          foundMe=true;
        for (int j=0;j<magDavid.size();j++){
          if ((fabs(MagK(i)-magDavid[j])<1e-5)){
            if (!foundMe){
              foundMe=true;
              for (int nL = NumImages; nL <= NumImages; nL++) {
                double kdotr;
                double rp;
                if (diagonal){
                  rp=(r+nL*box[0])/sqrt(3);
                  kdotr=(rp*kVecs(i)[0]+rp*kVecs(i)[1]+rp*kVecs(i)[2]);
                }
                else {
                  rp=r+nL*box[0];
                  kdotr=rp*kVecs(i)[0];
                }
                V+=cos(kdotr)*uk[j]/(box[0]*box[1]*box[2]); //*uk;
              }
            }
          }
        }
        if (!foundMe){
          cerr<<MagK[i]<<" "<<endl;
          cerr<<"ERROR! "<<endl;
          exit(1);
        }
      }
      outFile<<r<<" "<<V<<" "<<VShort_r+Vk0<<endl;
    }
  }

  void BasicEwald(double &alpha)
  {
    Array<dVec,1> kVecs;
    Array<double,1> MagK;
    SetupkVecs3D(kVecs,MagK,true);
    double vol=box[0]*box[1]*box[2];
    ofstream outfile;

    // k space part
    outfile.open("kData.txt");
    double Vlong_k0 = -4.0*M_PI*Z1Z2/(4.0*alpha*alpha*vol);
    cerr << "Vlong_k0 = " << Vlong_k0 << endl;
    outfile<<0.0<<" "<<Vlong_k0<<endl;
    vector<pair<double, double> > k_ewald;
    for (int i=0;i<kVecs.size();i++){
      dVec k=kVecs(i);
      double k2=dot(k,k);
      double kMag=sqrt(k2);
      double val=(4*M_PI*Z1Z2)/(vol*k2)*exp(-k2/(4*alpha*alpha));
      k_ewald.push_back(make_pair(kMag,val));
    }
    sort(k_ewald.begin(),k_ewald.end());
    vector<pair<double, double> >::iterator new_end = unique(k_ewald.begin(), k_ewald.end(), same);
    // delete all elements past new_end
    k_ewald.erase(new_end, k_ewald.end());
    for (int i=1;i<k_ewald.size();i++){
      outfile<<k_ewald[i].first<<" "<<k_ewald[i].second<<endl; // *vol<<endl;
    }
    outfile.close();

    // r space part
    outfile.open("rData.txt");
    double rmax = 0.75 * sqrt (dot(box,box));
    Array<double,1> Vshort_r(NumPoints), Vlong_r(NumPoints);
    LinearGrid LongGrid;
    LongGrid.Init (r0, rmax, NumPoints);
    for (int i=0;i<LongGrid.NumPoints;i++){
      double r=LongGrid(i);
      Vshort_r(i) = 0.0;
      for (int nL = -NumImages; nL <= NumImages; nL++) {
        double rp = r + nL*box[0];
        Vshort_r(i) += Z1Z2*erfc(alpha*rp)/rp;
        Vlong_r(i) += Z1Z2*erf(alpha*rp)/rp;
      }
      outfile<<r<<" "<<Vshort_r(i)<<" "<<Vlong_r(i)<<endl;
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
  e.SetkCut(atof(argv[2]));
  e.Setr0(atof(argv[3]));
  double rCut = atof(argv[4]);
  e.SetrCut(rCut);
  e.SetCharge(atof(argv[5]));
  e.SetNumPoints(atoi(argv[6]));
  bool doBreakup = atoi(argv[7]);

  if (doBreakup) {
    e.Breakup();
    e.SetNumImages(0);
  } else {
    double alpha=10.0/(L*2);
    alpha = rCut/L;
    e.SetNumImages(0);
    e.BasicEwald(alpha);
  }
  bool diagonal = atoi(argv[8]);
  e.FourierTransform(diagonal);
}
