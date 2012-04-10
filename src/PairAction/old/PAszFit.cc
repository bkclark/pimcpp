#include "PAFit.h"

#ifdef MAKE_FIT

void PAszFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  GridIsMine = true;
  assert(inSection.ReadVar ("Order", Order));
  UsePBC = inSection.ReadVar ("Box", Box);
}


void PAszFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "szfit");
  outSection.WriteVar ("Order", Order);
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
}


void PAszFitClass::DoFit (Rho &rho)
{
  int NumCoefs = (Order+1)*(Order+2)/2 - 1;
  NumBetas = 1;
  Array<double,2>  UCoefs(qgrid->NumPoints, NumCoefs);
  Array<double,2> dUCoefs(qgrid->NumPoints, NumCoefs);
  UCoefs = 0.0;
  dUCoefs = 0.0;

  Ukj.resize(1);
  dUkj.resize(1);
  sMax.resize(1);

  Ukj(0).Init(qgrid, UCoefs);
  dUkj(0).Init(qgrid, dUCoefs);

}

void PAszFitClass::WriteFit (IOSectionClass &outSection)
{
  int NumCoefs = (Order+1)*(Order+2)/2 - 1;
  Array<double,2> UCoefs(qgrid->NumPoints, NumCoefs); 
  Array<double,2> dUCoefs(qgrid->NumPoints, NumCoefs); 
  double beta = SmallestBeta;
  outSection.NewSection("Fit");
  outSection.WriteVar ("beta", beta);
  for (int j=0; j<qgrid->NumPoints; j++)
    for (int k=0; k<NumCoefs; k++) {
      UCoefs(j,k) = Ukj(i)(j,k);
      dUCoefs(j,k) = Ukj(i)(j,k);
    }
  outSection.WriteVar ("UCoefs", UCoefs);
  outSection.WriteVar ("dUCoefs", dUCoefs);
  outSection.WriteVar ("sMax", sMax(0));
  outSection.CloseSection();
}

void PAszFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight = 0.0;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *tdat = fopen ("t.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *ydat = fopen ("y.dat", "w");
  LinearGrid qgrid2(qgrid->Start, 0.999*qgrid->End, 20);
  Array<double,1> Ul, dUl;
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    // HACK
    //q = 1.0;
    double zmax = 0.9999999*min(2.0*q,sMax(level));
    LinearGrid zgrid(0.0, zmax, 20);
    for (int zi=0; zi<zgrid.NumPoints; zi++) {
      double z = zgrid(zi);
      double y = z/zmax;
      double smax = 0.9999999*min(2.0*q,sMax(level));
      //cerr << "smin = "  << z << " smax = " << smax << endl;
      LinearGrid sgrid(z, smax, 100);
      for (int si=0; si<sgrid.NumPoints; si++) {
	double s = sgrid(si);
	//cerr << "q = " << q << " z = " << z << " s = " << s << endl;
	double t = s/smax;
	double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
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

	rho.U_lArray(r,rp,Ul, dUl);
	rho.UdU(r, rp, costheta, Ul, dUl, Uex, dUex);
	if (!isnan(Uex) && !isnan(dUex)) {
	  Ufit = U(q, z, s*s, level);
	  dUfit = dU(q, z, s*s, level);
	  U2err += w*(Uex-Ufit)*(Uex-Ufit);
	  dU2err += w*(dUex-dUfit)*(dUex-dUfit);
	  weight += w;
	}
	fprintf (Uxdat, "%1.16e ", Uex);
	fprintf (Ufdat, "%1.16e ", Ufit);
	fprintf (dUxdat, "%1.16e ", dUex);
	fprintf (dUfdat, "%1.16e ", dUfit);
	fprintf (tdat, "%1.16e ", t);
	fprintf (costhetadat, "%1.16e ", costheta);
	fprintf (ydat, "%1.16e ", y);
      }
      fprintf (Uxdat, "\n");
      fprintf (Ufdat, "\n");
      fprintf (dUxdat, "\n");
      fprintf (dUfdat, "\n");
      fprintf (tdat, "\n");
      fprintf (ydat, "\n");
      fprintf (costhetadat, "\n");
    }
  }
  fclose (Uxdat); fclose(Ufdat); fclose(tdat); fclose(ydat); 
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}

#endif


bool PAszFitClass::Read (IOSectionClass &in,
			 double smallestBeta, int numBetas)
{
  NumBetas = numBetas;
  SmallestBeta = smallestBeta;
  // Resize
  //  Usplines.resize(NumBetas);
  //  dUsplines.resize(NumBetas);
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

  assert (in.ReadVar("Order", Order));
  assert (in.OpenSection("qGrid"));
  qgrid = ReadGrid(in);
  in.CloseSection();
  int NumCoefs = (Order+1)*(Order+2)/2 - 1;

  Array<double,2> UCoefs(qgrid->NumPoints, NumCoefs); 
  Array<double,2> dUCoefs(qgrid->NumPoints, NumCoefs); 

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
    assert (in.ReadVar ("UCoefs", UCoefs));
    assert (in.ReadVar ("dUCoefs", dUCoefs));
  }
  in.CloseSection(); // "Fits"

  return (true);
}

double PAszFitClass::U(double q, double z, double s2, int level)
{
  return (0.0);
}

double PAszFitClass::dU(double q, double z, double s2, int level)
{
  return (0.0);
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
//       double smax = min (2.0*q, sqrt(-4.0*lambda*myBeta*log(1.0e-4)));
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
