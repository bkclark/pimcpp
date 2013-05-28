#include "PAFit.h"
#include "../SpecialFunctions/HermitePoly.h"
#include "../SpecialFunctions/LegendrePoly.h"
#include "../SpecialFunctions/PolynomialSet.h"
#include "../Fitting/Fitting.h"

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAcoulombFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  GridIsMine = true;
  assert(inSection.ReadVar ("Order", Order));
  Ucoefs.resize(Order+1);
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAcoulombFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "coulombfit");
  outSection.WriteVar ("Order", Order);
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
}



// void PAcoulombFitClass::AddFit (Rho &rho)
// {
//   NumBetas++;
//   Uj.resizeAndPreserve(NumBetas);
//   dUj.resizeAndPreserve(NumBetas);

//   FILE *Udebug = fopen ("U.dat", "w");
//   FILE *sdebug = fopen ("s.dat", "w");

//   lambda = rho.lambda;
//   double beta = rho.Beta();
//   const int N = 400;
//   int M = Order+1;
//   Array<double,2>  UCoefs(qgrid->NumPoints, M);
//   Array<double,2> dUCoefs(qgrid->NumPoints, M);

//   Array<double,1> Uexact(N), dUexact(N), sigma(N);
//   Array<double,2> Basis(N,M-1);
//   Array<double,1> Hn(2*M), errors(M), Ui(M), dUi(M);
//   for (int qi=0; qi<qgrid->NumPoints; qi++) {
//     double q = (*qgrid)(qi);
//     LinearGrid sgrid(0.0, 2.0*q, N);
//     double U0, dU0;
//     rho.UdU_Coulomb (q, q, 1.0, U0, dU0);
//     for (int si=0; si<sgrid.NumPoints; si++) {
//       double s = sgrid(si);
//       double costheta;
//       if (q == 0.0)
// 	costheta = 1.0;
//       else
// 	costheta = 1.0 - s*s/(2.0*q*q);
//       costheta = min(1.0, costheta);
//       costheta = max(-1.0, costheta);
//       // Compute exact U and dU
//       double Uex, dUex;
//       rho.UdU_Coulomb (q, q, costheta, Uex, dUex);
//       Uexact(si) = Uex-U0;
//       fprintf (Udebug, "%1.16e ", Uex-U0);
//       fprintf (sdebug, "%1.16e ", s);
//       fflush (Udebug);
//       fflush (sdebug);
//       dUexact(si) = dUex-dU0;
//       if (isnan(Uexact(si))) {
// 	cerr << "NAN is Uexact for q = " << q << " s = " << s << endl;
// 	abort();
//       }
//       if (isnan(dUexact(si))){
// 	cerr << "NAN in dUexact for q = " << q << "s = " << s << endl;
// 	abort();
//       }
//       // Compute weight for point
//       double rho = exp(-s*s/(4.0*lambda*beta));
//       if (rho > 1e-8)
// 	sigma(si) = 1.0/sqrt(rho);
//       else
// 	sigma(si) = 1.0/sqrt(1e-8);
//       // Compute basis functions
//       double t = s / sqrt(4.0*lambda*beta);
//       for (int i=1; i<M; i++) 
// 	Basis(si,i-1) = pow(s, 2*i);
//     }
//     // Now do fits
//     LinFitSVD( Uexact, sigma, Basis,  Ui, errors, 1.0e-12);
//     LinFitSVD(dUexact, sigma, Basis, dUi, errors, 1.0e-12);
//     UCoefs(qi,0) = U0;
//     dUCoefs(qi,0) = dU0;
//     UCoefs(qi,Range(1,M-1)) = Ui;
//     dUCoefs(qi,Range(1,M-1)) = dUi;
//     fprintf (Udebug, "\n");
//     fflush (Udebug);
//     fprintf (sdebug, "\n");
//     fflush (sdebug);
//   } 
    
//   fclose (Udebug);
//   fclose (sdebug);

//   // Initialize splines
//   Uj(NumBetas-1).Init(qgrid, UCoefs);
//   dUj(NumBetas-1).Init(qgrid, dUCoefs);  
// }



class CoulombFitIntegrand
{
public:
  Rho &rho;
  int n;
  double q;
  double smax;

  double operator()(double x)
  {
    double s = 0.5*(x+1.0)*smax;

    double costheta;
    if (q == 0.0)
      costheta = 1.0;
    else
      costheta = 1.0 - s*s/(2.0*q*q);
    costheta = min(1.0, costheta);
    costheta = max(-1.0, costheta);

    double U, dU;
    rho.UdU_Coulomb(q,q,costheta, U, dU);
    return (0.5*(2.0*n+1)*U*LegendrePoly(n,x));
  }

  CoulombFitIntegrand(Rho &myrho) : rho(myrho)
  {  }
};




// void PAcoulombFitClass::AddFit (Rho &rho)
// {
//   NumBetas++;
//   Ujshort.resizeAndPreserve(NumBetas);
//   Ujlong.resizeAndPreserve(NumBetas);
//   dUjshort.resizeAndPreserve(NumBetas);
//   dUjlong.resizeAndPreserve(NumBetas);

//   lambda = rho.lambda;
//   double beta = rho.Beta();
//   int M = Order+1;
//   Array<double,2>  UCoefsShort(qgrid->NumPoints, M);
//   Array<double,2>  UCoefsLong(qgrid->NumPoints, M);
//   Array<double,2> dUCoefsShort(qgrid->NumPoints, M);
//   Array<double,2> dUCoefsLong(qgrid->NumPoints, M);

//   CoulombFitIntegrand integrand(rho);
//   const double Tolerance = 1.0e-7;
//   double sigma = sqrt(2.0*lambda*beta);
//   for (int qi=0; qi<qgrid->NumPoints; qi++) {
//     integrand.q = (*qgrid)(qi);
//     cerr << "qi = " << qi << " of " << qgrid->NumPoints << endl;
//     integrand.smax = min(2.0*integrand.q, 3.0*sigma);
//     for (int n=0; n<M; n++) {
//       integrand.n=n;
//       GKIntegration<CoulombFitIntegrand,GK15> Fit_gk(integrand);
//       Fit_gk.SetRelativeErrorMode();
//       UCoefsShort(qi,n) = Fit_gk.Integrate(-1.0, 1.0, Tolerance);
//     }
//     integrand.smax = 2.0*integrand.q;
//     for (int n=0; n<M; n++) {
//       integrand.n=n;
//       GKIntegration<CoulombFitIntegrand,GK15> Fit_gk(integrand);
//       Fit_gk.SetRelativeErrorMode();
//       UCoefsLong(qi,n) = Fit_gk.Integrate(-1.0, 1.0, Tolerance);
//     }
//   }

//   // Initialize splines
//   Ujshort(NumBetas-1).Init(qgrid, UCoefsShort);
//   dUjshort(NumBetas-1).Init(qgrid, dUCoefsShort);  
//   Ujlong(NumBetas-1).Init(qgrid, UCoefsLong);
//   dUjlong(NumBetas-1).Init(qgrid, dUCoefsLong);  
// }






class LConvertIntegrandClass
{
public:
  const PolynomialClass &P;
  int n;
  double smax;
  double operator()(double x)
  {
    double s = 0.5*(x+1.0)*smax;
    return (0.5*(2.0*n+1)*P(s)*LegendrePoly(n,x));
  }
  LConvertIntegrandClass(const PolynomialClass &p) : P(p)
  { }
};




Array<double,1> ConvertToLegendre(const PolynomialClass &P,
				  double smax)
{
  Array<double,1> coefs(P.Order()+1);
  coefs = 0.0;
  LConvertIntegrandClass integrand(P);
  integrand.smax = smax;
  for (int n=0; n<=P.Order(); n++) {
    integrand.n = n;
    GKIntegration<LConvertIntegrandClass, GK15> integrator(integrand);
    //integrator.SetRelativeErrorMode();
    coefs(n) = integrator.Integrate(-1.0, 1.0, 1.0e-7, 1.0e-7, false);
  }
  return coefs;
}


inline double factorial(double n)
{
  if (n<=1.0)
    return(1.0);
  else
    return (n*factorial(n-1.0));
}

class HConvertIntegrandClass
{
public:
  const PolynomialClass &P;
  int n;
  double Sqrt4LambdaBeta;
  double operator()(double x)
  {
    double s = x*Sqrt4LambdaBeta;
    double Ak = 1.0/((double)(1<<n)*factorial((double)n)*sqrt(M_PI));
    return (Ak*exp(-x*x)*P(s)*HermitePoly(n,x));
  }
  HConvertIntegrandClass(const PolynomialClass &p) : P(p)
  { }
};


Array<double,1> ConvertToHermite(const PolynomialClass &P,
				 double lambda, double beta)
{
  Array<double,1> coefs(P.Order()+1);
  coefs = 0.0;
  HConvertIntegrandClass integrand(P);
  integrand.Sqrt4LambdaBeta = sqrt(4.0*lambda*beta);
  for (int n=0; n<=P.Order(); n+=2) {
    integrand.n = n;
    GKIntegration<HConvertIntegrandClass, GK15> integrator(integrand);
    //integrator.SetRelativeErrorMode();
    coefs(n) = integrator.Integrate(-10.0, 10.0, 1.0e-8);
  }
  return coefs;
}

class GaussianWeightClass : public WeightFuncClass
{
public:
  double FourLambdaBetaInv, MinWeight;
  
  /// returns a Gaussian or the miniumum weight, whichever is greater
  double operator()(double s)
  {
    double gauss = exp(-s*s*FourLambdaBetaInv) + MinWeight;
    //return ( (gauss>MinWeight) ? gauss : MinWeight);
    return (gauss);
  }
  GaussianWeightClass (double lambda, double beta, double minWeight)
  {
    FourLambdaBetaInv = 1.0/(4.0*lambda*beta);
    MinWeight = minWeight;
  }
};


class UFitIntegrand
{
public:
  Rho &rho;
  PolynomialSetClass &Pset;
  WeightFuncClass &WF;
  int n;
  double q;

  double operator()(double s)
  {
    double costheta;
    if (q == 0.0)
      costheta = 1.0;
    else
      costheta = 1.0 - s*s/(2.0*q*q);
    costheta = min(1.0, costheta);
    costheta = max(-1.0, costheta);

    double U, dU;
    rho.UdU_Coulomb(q,q,costheta, U, dU);
    return (WF(s)*U*Pset(n,s));
  }

  UFitIntegrand(Rho &myrho, 
		PolynomialSetClass &mySet,
		WeightFuncClass &wf) : rho(myrho), Pset(mySet), WF(wf)
  {  }
};



void PAcoulombFitClass::AddFit (Rho &rho)
{
  NumBetas++;
  Ujshort.resizeAndPreserve(NumBetas);
  dUjshort.resizeAndPreserve(NumBetas);
  Ujlong.resizeAndPreserve(NumBetas);
  dUjlong.resizeAndPreserve(NumBetas);

  lambda = rho.lambda;
  double beta = rho.Beta();
  int M = Order+1;

  Array<double,2>  UCoefsShort(qgrid->NumPoints, M);
  Array<double,2> dUCoefsShort(qgrid->NumPoints, M);
  Array<double,2>  UCoefsLong(qgrid->NumPoints, M);
  Array<double,2> dUCoefsLong(qgrid->NumPoints, M);

  UCoefsShort = 0.0;  UCoefsLong=0.0;
  dUCoefsShort = 0.0; dUCoefsLong=0.0;

  const double Tolerance = 1.0e-8;
  GaussianWeightClass gaussWeight(lambda, beta, 1.0e-8);
  PolynomialClass poly(Order);
  PolynomialSetClass Pset;
  UFitIntegrand integrand(rho, Pset, gaussWeight);
  for (int qi=0; qi<qgrid->NumPoints; qi++) {
    cerr << "qi = " << qi+1 << " of " << qgrid->NumPoints << endl;
    integrand.q = (*qgrid)(qi);
    double smax = 2.0*integrand.q;
    // Create the orthonormal polynomials
    Pset.MakeOrthoSet(Order, -smax, smax, gaussWeight);
    for (int j=0; j<M; j++)
      poly[j] = 0.0;
    cerr << "Fitting to custom polynomial set:\n";
    for (int n=0; n<M; n+=2) {
      integrand.n=n;
      GKIntegration<UFitIntegrand,GK15> Fit_gk(integrand);
      Fit_gk.SetRelativeErrorMode();
      double C = Fit_gk.Integrate(-smax, smax, Tolerance);
      //poly[n] = C;
      //UCoefs(qi, n) = C;
      for (int j=0; j<M; j++)
 	poly[j] += C * Pset(n)[j];
    }
    cerr << "Hermite conversion:\n";
    UCoefsLong(qi, Range::all()) = ConvertToHermite(poly, lambda, beta);
    cerr << "Legendre conversion:\n";
    UCoefsShort(qi, Range::all()) = ConvertToLegendre(poly, 2.0*integrand.q);
  }

  // Initialize splines
  Ujshort(NumBetas-1).Init(qgrid, UCoefsShort);
  dUjshort(NumBetas-1).Init(qgrid, dUCoefsShort);  
  Ujlong(NumBetas-1).Init(qgrid, UCoefsLong);
  dUjlong(NumBetas-1).Init(qgrid, dUCoefsLong);  

}



void PAcoulombFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *sdat = fopen ("s.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *qdat = fopen ("q.dat", "w");
  LinearGrid qgrid2(qgrid->Start, qgrid->End, 315);
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    LinearGrid sgrid(0.0, 2.0*q, 350);
    for (int si=0; si<sgrid.NumPoints; si++) {
      double s = sgrid(si);
      double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
      double Uex, dUex, Ufit, dUfit;
      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q); 
      rho.UdU_Coulomb(q, q, costheta, Uex, dUex);
      Ufit = U(q, 0.0, s*s, level);
      dUfit = dU(q, 0.0, s*s, level);
      U2err += w*(Uex-Ufit)*(Uex-Ufit);
      dU2err += w*(dUex-dUfit)*(dUex-dUfit);
      weight += w;
      fprintf (Uxdat, "%1.16e ", Uex);
      fprintf (Ufdat, "%1.16e ", Ufit);
      fprintf (sdat, "%1.16e ", s);
      fprintf (costhetadat, "%1.16e ", costheta);
      fprintf (qdat, "%1.16e ", q);
    }
    fprintf (Uxdat, "\n");
    fprintf (Ufdat, "\n");
    fprintf (sdat, "\n");
    fprintf (qdat, "\n");
    fprintf (costhetadat, "\n");
  }
  fclose (Uxdat); fclose(Ufdat); fclose(sdat); fclose(qdat); 
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void PAcoulombFitClass::WriteFits (IOSectionClass &outSection)
{
  Array<double,2> UCoefsShort(qgrid->NumPoints, Order+1); 
  Array<double,2> dUCoefsShort(qgrid->NumPoints, Order+1); 
  Array<double,2> UCoefsLong(qgrid->NumPoints, Order+1); 
  Array<double,2> dUCoefsLong(qgrid->NumPoints, Order+1); 
  double beta = SmallestBeta;
  for (int i=0; i<NumBetas; i++) {
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int j=0; j<qgrid->NumPoints; j++)
      for (int k=0; k<(Order+1); k++) {
	UCoefsShort(j,k) = Ujshort(i)(j,k);
	dUCoefsShort(j,k) = dUjshort(i)(j,k);
	UCoefsLong(j,k) = Ujlong(i)(j,k);
	dUCoefsLong(j,k) = dUjlong(i)(j,k);
      }
    outSection.WriteVar ("Ujshort", UCoefsShort);
    outSection.WriteVar ("dUjshort", dUCoefsShort);
    outSection.WriteVar ("Ujlong", UCoefsLong);
    outSection.WriteVar ("dUjlong", dUCoefsLong);
    outSection.CloseSection();
    beta *= 2.0;
  }
}

#endif

// double PAcoulombFitClass::U(double q, double z, double s2, int level)
// {
//   if (q < qgrid->End) {
//     Uj(level)(q, Ucoefs);
//     double s2j = 1.0;
//     double Usum = 0.0;
//     for (int j=0; j<=Order; j++){
//       Usum += Ucoefs(j)*s2j;
//       s2j *= s2;
//     }
//     return (Usum);
//   }
//   else {
//     double beta = SmallestBeta;
//     for (int i=0; i<level; i++)
//       beta *= 2.0;
//     // Coulomb action is independent of z
//     return (beta*Pot->V(q));
//   }
// }

inline double SmoothStep(double x)
{
  double a = 25.0;
  return (1.0/(1.0+exp(a*(x-1.0))));
}

double PAcoulombFitClass::U(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;

  if (q < qgrid->End) {
    Array<double,1> Pn(Order+1);
    double s = sqrt(s2);
    double sigma = sqrt(2.0*lambda*beta);
    double smax = 2.0*q;//min(2.0*q, 3.0*sigma);
    double x = 2.0*s/smax - 1.0;
    LegendrePoly(x, Pn);
    Ujshort(level)(q, Ucoefs);
    double UsumShort = 0.0;
    for (int j=0; j<=Order; j++) 
      UsumShort += Ucoefs(j)*Pn(j);
    
    smax = 2.0*q;
    x = s*sqrt(4.0*lambda*beta);
    HermitePoly(x, Pn);
    double UsumLong = 0.0;
    Ujlong(level)(q, Ucoefs);
    for (int j=0; j<=Order; j++) 
      UsumLong += Ucoefs(j)*Pn(j);

    double y = q/(2.0*sigma);
    double f = SmoothStep(y);

    return (f*UsumShort+(1.0-f)*UsumLong);
  }
  else {
    // Coulomb action is independent of z
    return (beta*Pot->V(q));
  }
}


// double PAcoulombFitClass::U(double q, double z, double s2, int level)
// {
//   if (q < qgrid->End) {
//     double s = sqrt(s2);
//     Uj(level)(q, Ucoefs);
//     double Usum = 0.0;
//     double s2j = 1.0;
//     for (int j=0; j<=Order; j++) {
//       Usum += Ucoefs(j)*s2j;
//       s2j *= s;
//     } 
//     return (Usum);
//   }
//   else {
//     double beta = SmallestBeta;
//     for (int i=0; i<level; i++)
//       beta *= 2.0;
//     // Coulomb action is independent of z
//     return (beta*Pot->V(q));
//   }
// }

double PAcoulombFitClass::dU(double q, double z, double s2, int level)
{
  if (q < qgrid->End) {
    dUjshort(level)(q, Ucoefs);
    double s2j = 1.0;
    double dUsum = 0.0;
    for (int j=0; j<=Order; j++){
      dUsum += Ucoefs(j)*s2j;
      s2j *= s2;
    }
    return (dUsum);
  }
  else
    return 0.0;
}




bool PAcoulombFitClass::Read (IOSectionClass &in,
			      double smallestBeta, int NumBetas)
{
  // Resize
  Ujshort.resize(NumBetas);
  dUjshort.resize(NumBetas);
  Ujlong.resize(NumBetas);
  dUjlong.resize(NumBetas); 
  Array<double,2> temp;

  // Read Particles;
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;

  // Read Potential;
  assert(in.OpenSection("Potential"));
  Pot = ReadPotential(in);
  in.CloseSection();

  // Read the fits
  assert(in.OpenSection("Fits"));
  // Read the qgrid
  assert(in.OpenSection("qGrid"));
  qgrid = ReadGrid (in);
  in.CloseSection();
  GridIsMine=true;
  // Read Order
  assert (in.ReadVar ("Order", Order));
  Ucoefs.resize(Order+1);

  double desiredBeta = smallestBeta;
  for (int betaIndex=0; betaIndex<NumBetas; betaIndex++) {
    int i=0;
    int numTemps = in.CountSections("Fit");
    bool found=false;
    while ((i<numTemps) && !found) {
      double beta;
      assert (in.OpenSection("Fit", i));
      assert (in.ReadVar("beta", beta));
      if ((fabs(beta-desiredBeta)/desiredBeta) < 1.0e-12)
	found = true;
      else {
	in.CloseSection();
	i++;
      }
    }
    if (!found) {
    cerr << "Couldn't find beta = " << desiredBeta 
	 << " in fit file.  Exitting\n";
    exit(1);
    }
    // Now read the fit coefficents
    assert(in.ReadVar("Ujshort", temp));
    Ujshort(betaIndex).Init(qgrid,temp);
    assert(in.ReadVar("Ujlong", temp));
    Ujlong(betaIndex).Init(qgrid,temp);
    assert(in.ReadVar("dUjshort", temp));
    dUjshort(betaIndex).Init(qgrid,temp);
    assert(in.ReadVar("dUjlong", temp));
    dUjlong(betaIndex).Init(qgrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



// void Rho::WriteFit (string fileName)
// {
//   FILE *fout = fopen (fileName.c_str(), "w");
//   FILE *fpout = fopen ("var.dat", "w");
//   assert (fout != NULL);
//   int lmax = U_ls.size()-1;
//   LinearGrid qgrid (1e-4, 7.5, 101);
  
//   // HACK
//   double myBeta = 0.125;
//   //double myBeta = Beta();

//   for (int qindex=0; qindex<qgrid.NumPoints; qindex++)
//     {
//       cerr << " q = " << qindex << endl;
//       double q = qgrid(qindex);
//       double min (2.0*q, sqrt(-4.0*lambda*myBeta*log(1.0e-4)));
//       LinearGrid zgrid(-smax, smax, 101);
//       LinearGrid sgrid(0.0, smax, 101);
//       for (int zindex=0; zindex<zgrid.NumPoints; zindex++)
// 	{
// 	  double z = zgrid(zindex);
// 	  double r = q + 0.5*z;
// 	  double rp = q - 0.5*z;
// 	  if (r <= 0.0) r = 1.0e-5;
// 	  if (rp <= 0.0) rp = 1.0e-5;
// 	  double x = Transform.r2x(r);
// 	  double xp = Transform.r2x(rp);
// 	  Array<double,1> Ulvec(lmax+1), dUlvec(lmax+1);
// 	  // First, calculate the diagonal parts
// 	  U_lArray(r,r, Ulvec, dUlvec);
// 	  double Ur, Urp, dUr, dUrp;
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// 	  U_lArray(rp,rp, Ulvec, dUlvec);
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,x);
// // 	  double Udiag_r = U(r,r,1.0, Ularray);

// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(xp,xp);
// // 	  double Udiag_rp = U(rp,rp,1.0, Ularray);

// 	  // Now calculate off-diagonal elements
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,xp);
// 	  U_lArray(r,rp, Ulvec, dUlvec);
// 	  //LinearGrid sgrid(fabs(r-rp)+1e-9, r+rp-1e-9, 101);
// 	  for (int sindex=0; sindex<sgrid.NumPoints; sindex++)
// 	    {
// 	      double s = sgrid(sindex);	      
// 	      double costheta;
// 	      if ((r<=0.0) || (rp<=0.0))
// 		costheta = 1.0;
// 	      else
// 		costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
// 	      //if (costheta > 1.0)  costheta = 1.0;
// 	      //if (costheta < -1.0) costheta = -1.0;
// 	      double theta = acos (costheta);
// 	      double Uval, dUval, Ufitval, dUfitval, Ufit2val, dUfit2val;
// 	      if (isnan(theta))
// 		{
// // 		  cerr << "NAN in theta.\n"
// // 		       << "costheta = " << costheta << endl
// // 		       << "r = " << r << " rp = " << rp << endl 
// // 		       << "q = " << q << " z = " << z << " s = " << s << endl;
// 		  Uval = NAN;	   dUval = NAN;
// 		  Ufitval = NAN;   dUfitval = NAN;
// 		  Ufit2val = NAN;  dUfit2val = NAN;
// 		}
// 	      else
// 		{
// 		  if ((r < grid->End) && (rp < grid->End)) {
// 		    UdU(r,rp,costheta, Ulvec, dUlvec, Uval, dUval);
// 		    //Uval = U(r,rp,costheta, Ularray);
// 		    //Uval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufitval = Ufit(q, s, z);
// 		    //Ufitval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufit2val = Ufit2(q,s,z);
// 		    dUfit2val = Ufit2.dU(q,s,z);
// 		    //Ufit2val -= 0.5*(Udiag_r + Udiag_rp);
// 		  }
// 		  else {
// 		    cerr << "r or rp outside grid.\n";
// 		    Uval = 0.0;
// 		  }
// 		}
// 	      fprintf (fout, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e\n", 
// 		       Uval, Ufitval, Ufit2val, dUval, dUfitval, dUfit2val);
// 	      fprintf (fpout, "%1.6e %1.6e %1.6e %1.8e\n", q, z, s,
// 		       exp(-s*s/(2.0*beta)));
// 	    }
// 	}
//     }
//   fclose (fout);
//   fclose(fpout);
// }
