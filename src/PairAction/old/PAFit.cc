#include "PAFit.h"

void PAszFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "szfit");
  outSection.WriteVar ("Order", Order);
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
}


void PAszFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  GridIsMine = true;
  assert(inSection.ReadVar ("Order", Order));
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAszFitClass::AddFit (Rho &rho)
{
  int NumCoefs = (Order+1)*(Order+2)/2 - 1;
  NumBetas++;
  Array<double,2>  UCoefs(qgrid->NumPoints, NumCoefs);
  Array<double,2> dUCoefs(qgrid->NumPoints, NumCoefs);
  UCoefs = 0.0;
  dUCoefs = 0.0;

  Ukj.resizeAndPreserve(NumBetas);
  dUkj.resizeAndPreserve(NumBetas);
  smax.resizeAndPreserve(NumBetas);

  Ukj(NumBetas-1).Init(qgrid, UCoefs);
  dUkj(NumBetas-1).Init(qgrid, dUCoefs);

}

void PAszFitClass::WriteFits (IOSectionClass &outSection)
{
  int NumCoefs = (Order+1)*(Order+2)/2 - 1;
  Array<double,2> UCoefs(qgrid->NumPoints, NumCoefs); 
  Array<double,2> dUCoefs(qgrid->NumPoints, NumCoefs); 
  double beta = SmallestBeta;
  for (int i=0; i<NumBetas; i++) {
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int j=0; j<qgrid->NumPoints; j++)
      for (int k=0; k<NumCoefs; k++) {
	UCoefs(j,k) = Ukj(i)(j,k);
	dUCoefs(j,k) = Ukj(i)(j,k);
      }
    outSection.WriteVar ("Ucoefs", UCoefs);
    outSection.WriteVar ("dUcoefs", dUCoefs);
    outSection.WriteVar ("smax", smax(i));
    outSection.CloseSection();
    beta *= 2.0;
  }
}

bool PAszFitClass::Read (IOSectionClass &inSection,
			 double lowestBeta, int NumBetas)
{

}

double PAszFitClass::U(double r, double rp, double costheta, int level)
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
