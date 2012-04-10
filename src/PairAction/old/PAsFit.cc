#include "PAFit.h"
#include "../Splines/BicubicSpline.h"
#include "../SpecialFunctions/LegendrePoly.h"
#include "../SpecialFunctions/PolynomialSet.h"
#include "../Distributed/DistributedMat.h"


const double Rho0Min = 1.0e-3;

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAsFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("yGrid"));
  ygrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.ReadVar("Order", Order));
  Coefs.resize(Order+1);
  Pn.resize(Order+1);
  GridsAreMine = true;
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAsFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    outSection.WriteVar ("Type", "sfit");
    outSection.NewSection("qGrid");
    qgrid->Write(outSection);
    outSection.CloseSection();
    outSection.NewSection("yGrid");
    ygrid->Write(outSection);
    outSection.CloseSection();
    outSection.WriteVar("Order", Order);
  }
}


// class sFitIntegrand
// {
// private:
//   double q, z;
//   double r, rp;
//   double sMin, sMax;
// public:
//   Rho &rho;
//   bool IsdU;
//   int k;

//   Array<double,1> Ul, dUl;
  
//   void Setqz(double newq, double newz,
// 	     double smin, double smax)
//   {
//     q = newq; z = newz;
//     sMin = smin; sMax = smax;
    
//     r = q + 0.5*z;
//     rp = q - 0.5*z;
//     rho.U_lArray(r,rp, Ul, dUl);
//   }

//   inline double operator()(double x)
//   {
//     double s=sMin+0.5*(sMax-sMin)*(x+1.0);
//     double costheta;
    
//     if ((r*rp)==0.0)
//       costheta = 1.0;
//     else
//       costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
    
//     costheta = min(1.0, costheta);
//     costheta = max(-1.0, costheta);
//     double U, dU;
//     rho.UdU(r,rp,costheta, Ul, dUl, U, dU);
//     if (isnan(U)) {
//       cerr << "r = " << r << " rp = " << rp << endl;
//       cerr << "costheta = " << costheta << endl;
//       cerr << "x = " << x << endl;
//     }
//     // U = 0.0;
//     //cerr << "r = " << r << " rp = " << rp << endl;
//     //cerr << "smin = " << sMin << " smax = " << sMax << endl;
//     //cerr << "x = " << x << endl;
//     //cerr << "costheta = " << costheta << endl;
//     if (!IsdU)
//       return (0.5*(2.0*k+1.0)*U*LegendrePoly(k,x));
//     else
//       return (0.5*(2.0*k+1.0)*dU*LegendrePoly(k,x));
//   }
  
//   sFitIntegrand(Rho &myRho) : rho(myRho)
//   {
//     Ul.resize(rho.U_ls.size());
//     dUl.resize(rho.U_ls.size());
//   };
// };



class sFitIntegrand
{
private:
  double q, z;
  double r, rp;
  double sMin, sMax;
public:
  Rho &rho;
  PolynomialSetClass &Pset;
  WeightFuncClass &Weight;
  bool IsdU;
  int k;

  Array<double,1> Ul, dUl;
  
  void Setqz(double newq, double newz,
	     double smin, double smax)
  {
    q = newq; z = newz;
    sMin = smin; sMax = smax;
    
    r = q + 0.5*z;
    rp = q - 0.5*z;
    rho.U_lArray(r,rp, Ul, dUl);
  }

  inline double operator()(double s)
  {
    double costheta;
    
    if ((r*rp)==0.0)
      costheta = 1.0;
    else
      costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
    
    costheta = min(1.0, costheta);
    costheta = max(-1.0, costheta);
    double U, dU;
    rho.UdU(r,rp,costheta, Ul, dUl, U, dU);
    if (isnan(U) || isnan(dU)) {
      cerr << "r = " << r << " rp = " << rp << endl;
      cerr << "costheta = " << costheta << endl;
      cerr << "s = " << s << endl;
    }
    double p = Pset(k,s);
    double w = Weight(s);
    if (isnan(p))
      cerr << "NAN in p.\n";
    if (isnan(w))
      cerr << "NAN in w.\n";
    if (!IsdU)
      return (Weight(s)*Pset(k,s)*U);
    else
      return (Weight(s)*Pset(k,s)*dU);
  }
  
  sFitIntegrand(Rho &myRho, PolynomialSetClass &pset,
		WeightFuncClass &weight) : rho(myRho), Pset(pset),
					   Weight(weight)
  {
    Ul.resize(rho.U_ls.size());
    dUl.resize(rho.U_ls.size());
  }
};


class GaussianClass : public WeightFuncClass
{
public:
  double FourLambdaBetaInv;
  double operator()(double s)
  {
    return (exp(-s*s*FourLambdaBetaInv));
  }
  GaussianClass (double lambda, double beta)
  {
    FourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  }

};



/// y = |z/z_max|
/// z_max = 
void PAsFitClass::AddFit (Rho &rho)
{
  const double Tolerance = 1.0e-6;
  NumBetas++;
  Usplines.resizeAndPreserve(NumBetas);
  dUsplines.resizeAndPreserve(NumBetas);
  sMax.resizeAndPreserve(NumBetas);

  double lambda = Particle1.lambda + Particle2.lambda;
  double beta = rho.Beta();
  sMax(NumBetas-1) = sqrt(-4.0*lambda*beta*log(Rho0Min));

  GaussianClass weight(lambda,beta);
  PolynomialSetClass pset;

  int numq = qgrid->NumPoints;
  int numy = ygrid->NumPoints;
  CommunicatorClass comm;
  comm.SetWorld();
  DistributedArray3 Umat(numq, numy,Order+1,comm), 
    dUmat(numq, numy, Order+1, comm);
  sFitIntegrand integrand(rho, pset, weight);
  int qi, yi;
  qi = -1;
  for (int i=0; i<Umat.MyNumElements(); i++) {
    int tempqi;
    Umat.MyElement(i, tempqi, yi);
    if (tempqi != qi) {
      qi = tempqi;
      cerr << "qi = " << (qi+1) << " of " << numq << ".\n";      
    }
    cerr << "yi = " << yi << endl;
    double q = (*qgrid)(qi);    
    double zmax = 0.9999*min(2.0*q, sMax(NumBetas-1));
    double z = zmax*(*ygrid)(yi);
    double smin = z;
    double smax = 0.99999*min (2.0*q, sMax(NumBetas-1));
    //cerr << "<MakeOrthSet>\n";
    pset.MakeOrthoSet(Order, smin, smax, weight);
    //cerr << "</MakeOrthSet>\n";
    //cerr << "smin = " << smin << " smax = " << smax << endl;
    integrand.Setqz(q,z,smin,smax);
    PolynomialClass Upoly(Order), dUpoly(Order);
    Upoly = 0.0;
    dUpoly = 0.0;
    for(int k=0; k<=Order; k++) {
      integrand.k=k;
      integrand.IsdU=false;

      GKIntegration<sFitIntegrand,GK15> Uintegrator(integrand);
      // Accept a relative OR absolute tolerance of 1.0e-7
      double Uk = Uintegrator.Integrate(smin, smax,
					beta*Tolerance, Tolerance, 
					false);
      Upoly += (Uk*integrand.Pset(k));

      integrand.IsdU = true;
      // begin HACK
//       if (yi == 29)
// 	{
// 	  FILE *fout;
// 	  LinearGrid sg(smin, smax, 200);
// 	  fout = fopen ("stest.dat", "w");
// 	  for (int si=0; si<sg.NumPoints; si++)
// 	    fprintf (fout, "%1.12e %1.16e\n", sg(si), integrand(sg(si)));
// 	  fclose (fout);
// 	}
      // End HACK

      GKIntegration<sFitIntegrand,GK15> dUintegrator(integrand);
      double dUk = dUintegrator.Integrate(smin, smax, Tolerance,
					  Tolerance, false);
      dUpoly += (dUk*integrand.Pset(k));
    }
    for (int k=0; k<=Order; k++) {
      Umat(qi,yi,k) = Upoly[k];
      dUmat(qi,yi,k) = dUpoly[k];
    }
  }
  Umat.AllGather();
  dUmat.AllGather();
  // Initialize the bicubic splines
  Usplines(NumBetas-1).Init(qgrid,ygrid,Umat.Mat);
  dUsplines(NumBetas-1).Init(qgrid,ygrid,dUmat.Mat);
}



void PAsFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight = 0.0;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *sdat = fopen ("s.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *ydat = fopen ("y.dat", "w");
  LinearGrid qgrid2(qgrid->Start, qgrid->End, 21);
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    // HACK
    //q = 1.0;
    double zmax = 0.9999*min(2.0*q,sMax(level));
    LinearGrid zgrid(0.0, zmax, 20);
    for (int zi=0; zi<zgrid.NumPoints; zi++) {
      double z = zgrid(zi);
      double y = z/zmax;
      double smax = 2.0*q;//min(2.0*q,sMax(level));
      //cerr << "smin = "  << z << " smax = " << smax << endl;
      LinearGrid sgrid(z, smax, 100);
      for (int si=0; si<sgrid.NumPoints; si++) {
	double s = sgrid(si);
	double lambda = Particle1.lambda + Particle2.lambda;
	double w = exp(-s*s/(4.0*lambda*rho.Beta()));
	double Uex, dUex, Ufit, dUfit;
	double r, rp, costheta;
	r  = q+0.5*z;
	rp = q-0.5*z;
	if (q == 0.0)
	  costheta = 1.0;
	else
	  costheta = (r*r + rp*rp - s*s)/(2.0*r*rp); 
	
	//cerr << "costheta = " << costheta << endl;

	costheta = min(costheta,1.0);
	costheta = max(costheta,-1.0);


	rho.UdU(r, rp, costheta, Uex, dUex);
	Ufit = U(q, z, s*s, level);
	dUfit = dU(q, z, s*s, level);
	if (!isnan(Uex) && !isnan(dUex)) {
	  if (s <= sMax(level)) {
	    U2err += w*(Uex-Ufit)*(Uex-Ufit);
	    dU2err += w*(dUex-dUfit)*(dUex-dUfit);
	    weight += w;
	  }
	}
	fprintf (Uxdat, "%1.16e ", Uex);
	fprintf (Ufdat, "%1.16e ", Ufit);
	fprintf (dUxdat, "%1.16e ", dUex);
	fprintf (dUfdat, "%1.16e ", dUfit);
	fprintf (sdat, "%1.16e ", s);
	fprintf (costhetadat, "%1.16e ", costheta);
	fprintf (ydat, "%1.16e ", y);
      }
      fprintf (Uxdat, "\n");
      fprintf (Ufdat, "\n");
      fprintf (dUxdat, "\n");
      fprintf (dUfdat, "\n");
      fprintf (sdat, "\n");
      fprintf (ydat, "\n");
      fprintf (costhetadat, "\n");
    }
  }
  fclose (Uxdat); fclose(Ufdat); fclose(sdat); fclose(ydat); 
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void PAsFitClass::WriteFits (IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    Array<double,3> Umat(qgrid->NumPoints, ygrid->NumPoints, Order+1); 
    Array<double,3> dUmat(qgrid->NumPoints, ygrid->NumPoints,Order+1); 
    double beta = SmallestBeta;
    for (int bi=0; bi<NumBetas; bi++) {
      cerr << "Writing Beta "<< bi+1 << " of " << NumBetas << ":\n";
      outSection.NewSection("Fit");
      outSection.WriteVar ("beta", beta);
      double smax = sMax(bi);
      outSection.WriteVar ("sMax", smax);
      for (int qi=0; qi<qgrid->NumPoints; qi++)
	for (int yi=0; yi<ygrid->NumPoints; yi++) 
	  for (int j=0; j<=Order; j++) {
	    Umat(qi,yi,j)  =  Usplines(bi)(qi,yi,j);
	    dUmat(qi,yi,j) = dUsplines(bi)(qi,yi,j);
	  }
      outSection.WriteVar ("Umat", Umat);
      outSection.WriteVar ("dUmat", dUmat);
      outSection.CloseSection();
      beta *= 2.0;
    }
  }
}

#endif


// double PAsFitClass::U(double q, double z, double s2, int level)
// {
//   z = fabs(z);
//   double qmax = qgrid->End*1.000001;
//   double zmax = 1.000001*min (2.0*q, sMax(level));
//   double smax = 1.000001*min (2.0*q, sMax(level));
//   double smin = z;
//   double s=sqrt(s2);
//   double x;

//   if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
//     if (q == 0.0) {
//       Usplines(level)(0.0,0.0,Coefs);
//       x = 0.0;
//     }
//     else {
//       double y = z/zmax;
//       Usplines(level)(q,y, Coefs);
//       x = (s-smin)/(smax-smin)*2.0 - 1.0;
//     }

//     LegendrePoly(x, Pn);
//     // Now do summation
//     double sum=0.0;
//     for (int k=0; k<=Order; k++)
//       sum += Pn(k)*Coefs(k);
//     return (sum);
//   }
//   else {
//     double beta = SmallestBeta;
//     for (int i=0; i<level; i++)
//       beta *= 2.0;
//     double r = q+0.5*z;
//     double rp = q-0.5*z;
//     return (0.5*beta*(Pot->V(r)+Pot->V(rp)));
//   }
// }




double PAsFitClass::U(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End*1.00000001;
  double zmax = 0.9999*min (2.0*q, sMax(level));
  double s=sqrt(s2);
  double x;  

  if ((q<=qmax)&&(z<=zmax)) {
    if (q == 0.0) {
      Usplines(level)(0.0,0.0,Coefs);
    }
    else {
      double y = z/zmax;
      Usplines(level)(q,y, Coefs);
    }

    double s2n = 1.0;
    // Now do summation
    double sum=0.0;
    for (int k=0; k<=Order; k++) {
      sum += s2n*Coefs(k);
      s2n *= s;
    }
    return (sum);
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*beta*(Pot->V(r)+Pot->V(rp)));
  }
}


double PAsFitClass::dU(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End*1.00000001;
  double zmax = 0.9999*min (2.0*q, sMax(level));
  double s=sqrt(s2);
  double x;  

  if ((q<=qmax)&&(z<=zmax)) {
    if (q == 0.0) {
      dUsplines(level)(0.0,0.0,Coefs);
    }
    else {
      double y = z/zmax;
      dUsplines(level)(q,y, Coefs);
    }

    double s2n = 1.0;
    // Now do summation
    double sum=0.0;
    for (int k=0; k<=Order; k++) {
      sum += s2n*Coefs(k);
      s2n *= s;
    }
    return (sum);
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*(Pot->V(r)+Pot->V(rp)));
  }
}


// double PAsFitClass::dU(double q, double z, double s2, int level)
// {
//   z = fabs(z);
//   double qmax = qgrid->End*1.000001;
//   double zmax = 1.000001*min (2.0*q, sMax(level));
//   double smax = 1.000001*min (2.0*q, sMax(level));
//   double smin = z;
//   double s=sqrt(s2);
//   double x;

//   if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
//     if (q == 0) {
//       dUsplines(level)(0.0,0.0,Coefs);
//       x = 0.0;
//     }
//     else {
//       double y = z/zmax;
//       dUsplines(level)(q,y, Coefs);
//       x = (s-smin)/(smax-smin)*2.0 - 1.0;
//     }

//     LegendrePoly(x, Pn);
//     // Now do summation
//     double sum=0.0;
//     for (int k=0; k<=Order; k++)
//       sum += Pn(k)*Coefs(k);
//     return (sum);
//   }
//   else {
//     double beta = SmallestBeta;
//     for (int i=0; i<level; i++)
//       beta *= 2.0;
//     double r = q+0.5*z;
//     double rp = q-0.5*z;
//     return (0.5*(Pot->V(r)+Pot->V(rp)));
//   }
// }



bool PAsFitClass::Read (IOSectionClass &in,
			double smallestBeta, int numBetas)
{
  NumBetas = numBetas;
  SmallestBeta = smallestBeta;
  // Resize
  Usplines.resize(NumBetas);
  dUsplines.resize(NumBetas);
  Array<double,3> temp;

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
  assert(in.OpenSection("yGrid"));
  ygrid = ReadGrid (in);
  in.CloseSection();
  GridsAreMine=true;
  // Read Order
  assert(in.ReadVar("Order", Order));
  Coefs.resize(Order+1);
  Pn.resize(Order+1);
  sMax.resize(NumBetas);

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
    double smax;
    assert(in.ReadVar("sMax", smax));
    sMax(betaIndex) = smax;
    assert(in.ReadVar("Umat", temp));
    Usplines(betaIndex).Init(qgrid,ygrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,ygrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



