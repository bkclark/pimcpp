/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
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
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "DavidPAClass.h"

void DavidPAClass::ReadParams(IOSectionClass &in)
{

}

void DavidPAClass::WriteBetaIndependentInfo (IOSectionClass &out)
{

}

// void DavidPAClass::Error (Rho &rho, double &Uerror, double &dUerror)
// {

// }

// void DavidPAClass::DoFit (Rho &rho)
// {

// }

void DavidPAClass::WriteFit(IOSectionClass &out)
{


}

double DavidPAClass::dUdRTimesSigma(double r,int level)
{
  if (r>dUdRTimesSigmaSpline(level).grid->End){
    r=dUdRTimesSigmaSpline(level).grid->End;
  }
  if (r<dUdRTimesSigmaSpline(level).grid->Start){
    r=dUdRTimesSigmaSpline(level).grid->Start;
  }
  return dUdRTimesSigmaSpline(level)(r);
}


double DavidPAClass::dUdRTimesSigma_movers(double r,int level)
{
  if (r>dUdRTimesSigmaSpline_movers(level).grid->End){
    r=dUdRTimesSigmaSpline_movers(level).grid->End;
  }
  if (r<dUdRTimesSigmaSpline_movers(level).grid->Start){
    r=dUdRTimesSigmaSpline_movers(level).grid->Start;
  }
  return dUdRTimesSigmaSpline_movers(level)(r);
}


double DavidPAClass::d2UdR2TimesSigma(double r,int level)
{
  if (r>d2UdR2TimesSigmaSpline(level).grid->End){
    r=d2UdR2TimesSigmaSpline(level).grid->End;
  }
  if (r<d2UdR2TimesSigmaSpline(level).grid->Start){
    r=d2UdR2TimesSigmaSpline(level).grid->Start;
  }
  return d2UdR2TimesSigmaSpline(level)(r);
}

double DavidPAClass::d2UdR2TimesSigma_movers(double r,int level)
{
  if (r>d2UdR2TimesSigmaSpline_movers(level).grid->End){
    r=d2UdR2TimesSigmaSpline_movers(level).grid->End;
  }
  if (r<d2UdR2TimesSigmaSpline_movers(level).grid->Start){
    r=d2UdR2TimesSigmaSpline_movers(level).grid->Start;
  }
  return d2UdR2TimesSigmaSpline_movers(level)(r);
}


double DavidPAClass::U (double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  //  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  calcUsqzFast(s,q,z,level,uTemp);
  return uTemp;

}

double DavidPAClass::V(double r)
{
  if (r>ukj(0).grid->End){
    r=ukj(0).grid->End;
  }
  if (r<ukj(0).grid->Start){
    r=ukj(0).grid->Start;
  }
  return ukj(0)(0,r);
}

double DavidPAClass::dU(double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  //  cerr<<"my vtemp is "<<vTemp<<endl;
  return duTemp;


}
bool DavidPAClass::IsLongRange()
{
  return false;
}


void
DavidPAClass::Derivs (double q, double z, double s2, int level,
		      double &d_dq, double &d_dz)
{
  level += TauPos;
  
  double rmin = ukj(level).grid->Start;
  
  d_dq = 0.0;
  d_dq = 0.0;

  double r  = q + 0.5*z;
  double rp = q - 0.5*z;
  if (r  > ukj(level).grid->End)
    r  = ukj(level).grid->End;
  if (rp > ukj(level).grid->End)
    rp = ukj(level).grid->End;

  /////////////////////
  /// Diagonal part ///
  /////////////////////
  double dU_dr  = ukj(level).Deriv(1,r);
  double dU_drp = ukj(level).Deriv(1,rp);

  d_dq = 0.5*(dU_dr + dU_drp);
  d_dz = 0.25*(dU_dr - dU_drp);

  /////////////////////////
  /// Off-diagonal part ///
  /////////////////////////
  double z2=z*z;
  double s2inverse=1.0/s2;
  double Sto2k=s2;
  
  if (s2 > 0.0 && q<ukj(level).grid->End) { 
    (ukj(level))(q,TempukjArray); 
    (ukj(level)).Deriv(q,TempdukjArray); 
    for (int k=1; k<=NumTerms; k++) {  
      double Zto2j=1.0;
      double Zto2jm1=z;
      double currS=Sto2k;
      for (int j=0; j<=k; j++) {
	  // indexing into the 2darray
	double dq_coef  = TempdukjArray(k*(k+1)/2+j+1);
	double dz_coef  = TempukjArray (k*(k+1)/2+j+1);
	d_dq += dq_coef*Zto2j*currS;
	if (j > 0) {
	  d_dz += dz_coef*2.0*(double)j*Zto2jm1*currS;
	  Zto2jm1 *= z2;
	}
	Zto2j *= z2;
	currS=currS*s2inverse;				
      }				
      Sto2k=Sto2k*s2;
    } 
  }
}

void
DavidPAClass::DerivsFD (double q, double z, double s2, int level,
			 double &d_dq, double &d_dz)
{
  double s = sqrt(s2);
  double epsilon = 1.0e-6;

  double dummy;

  double qplus, qminus, zplus, zminus;
  calcUsqz(s, q+epsilon, z, level, qplus,  dummy, dummy);
  calcUsqz(s, q-epsilon, z, level, qminus, dummy, dummy);
  calcUsqz(s, q, z+epsilon, level, zplus,  dummy, dummy);
  calcUsqz(s, q, z-epsilon, level, zminus, dummy, dummy);
  d_dq = (qplus - qminus)/(2.0*epsilon);
  d_dz = (zplus - zminus)/(2.0*epsilon);
}

void
DavidPAClass::Derivs (double q, double z, double s2, int level,
		      double &d_dq, double &d_dz, double &d_ds)
{
  level += TauPos;
  
  double rmin = ukj(level).grid->Start;

  double r  = q + 0.5*z;
  double rp = q - 0.5*z;
  if (r  > ukj(level).grid->End)
    r  = ukj(level).grid->End;
  if (rp > ukj(level).grid->End)
    rp = ukj(level).grid->End;

  /////////////////////
  /// Diagonal part ///
  /////////////////////
  double dU_dr  = ukj(level).Deriv(1,r);
  double dU_drp = ukj(level).Deriv(1,rp);

  d_dq = 0.5*(dU_dr + dU_drp);
  d_dz = 0.25*(dU_dr - dU_drp);
  d_ds = 0.0;

  /////////////////////////
  /// Off-diagonal part ///
  /////////////////////////
  double z2=z*z;
  double s2inverse=1.0/s2;
  double sinverse = sqrt(s2inverse);
  double Sto2k=s2;
  
  if (s2 > 0.0 && q<ukj(level).grid->End) { 
    (ukj(level))(q,TempukjArray); 
    (ukj(level)).Deriv(q,TempdukjArray); 
    for (int k=1; k<=NumTerms; k++) {  
      double Zto2j=1.0;
      double Zto2jm1=z;
      double currS=Sto2k;
      for (int j=0; j<=k; j++) {
	  // indexing into the 2darray
	double dq_coef  = TempdukjArray(k*(k+1)/2+j+1);
	double dz_coef  = TempukjArray (k*(k+1)/2+j+1);
	d_dq += dq_coef*Zto2j*currS;
	d_ds += dz_coef*2.0*(double)(k-j)*Zto2j*currS*sinverse;
	if (j > 0) {
	  d_dz += dz_coef*2.0*(double)j*Zto2jm1*currS;
	  Zto2jm1 *= z2;
	}
	Zto2j *= z2;
	currS=currS*s2inverse;				
      }				
      Sto2k=Sto2k*s2;
    } 
  }
}

void
DavidPAClass::DerivsFD (double q, double z, double s2, int level,
			double &d_dq, double &d_dz, double &d_ds)
{
  double s = sqrt(s2);
  double epsilon = 1.0e-6;

  double dummy;

  double qplus, qminus, zplus, zminus, splus, sminus;
  calcUsqz(s, q+epsilon, z, level, qplus,  dummy, dummy);
  calcUsqz(s, q-epsilon, z, level, qminus, dummy, dummy);
  calcUsqz(s, q, z+epsilon, level, zplus,  dummy, dummy);
  calcUsqz(s, q, z-epsilon, level, zminus, dummy, dummy);
  calcUsqz(s+epsilon, q, z, level, splus,  dummy, dummy);
  calcUsqz(s-epsilon, q, z, level, sminus, dummy, dummy);
  d_dq = (qplus - qminus)/(2.0*epsilon);
  d_dz = (zplus - zminus)/(2.0*epsilon);
  d_ds = (splus - sminus)/(2.0*epsilon);
}

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void DavidPAClass::calcUsqz(double s,double q,double z,int level,
			    double &U, double &dU, double &V)
{
  //level=level+(NumTau-(TauPos+1));
  level=level+TauPos;
  //  cerr<<"My level is "<<level<<endl;
  
  double rmin = ukj(level).grid->Start;
  
  U=0.0;
  dU=0.0;
  //  level=level+4;
  // Check to make sure we're inside the grid.
  //  //  if (q > ukj(level).grid->End) {
  //  //    U = 0.0; dU=0.0; V = 0.0;
  //  //    return;
  //  //  }
  
  double r=q+0.5*z;
  double rprime=q-0.5*z;

 //  ///HACK! HACK! HACK! HACK!
//   if (r > ukj(level).grid->End) {
//     U =0.0 ; dU=0.0; V = 0.0;
//     return;
//   }
//     ///END HACK! END HACK! END HACK!

  if (r > ukj(level).grid->End) {
    //    U =0.0 ; dU=0.0; V = 0.0;
    //    return;
    r=ukj(level).grid->End;

  }
  if (rprime > ukj(level).grid->End) {
    //    U = 0.0; dU=0.0; V = 0.0;
    //    return;
    rprime=ukj(level).grid->End;

  }
//   if (q>ukj(level).grid->End){
//     q=ukj(level).grid->End;
//   }
//   q=0.5*(r+rprime);
//   z=z*0.5;
//   s=s*0.5;
 //   if (fabs(z)>q){
//     z=q;
//   }
//   if (s>q){
//     s=q;
//   }
  
  // This is the endpoint action 
  
  
  if (rprime < rmin) {
    //cerr << "rprime < rmin" << endl;
    rprime = rmin;
  }
  if(r < rmin) {
    //cerr << "r < rmin" << endl;
    r = rmin;
  }
  //if ((rprime < rmin) || (r < rmin)){
  //  //    cerr<<"I'm less then the min. Maybe this is messing me up\n";
  //  U = 5000.0; dU = 0.0;
  //  return;
  //}
  
//   if (q<ukj(level).grid->Start){
//     q=ukj(level).grid->Start;
//   }

  // Compensate for potential, which is subtracted from diaganal action in
  // dm file.
  V = 0.5*(ukj(level)(0,r) + ukj(level)(0,rprime));
  U+= 0.5*(ukj(level)(1,r)+ukj(level)(1,rprime));
  dU+=0.5*((dukj(level))(1,r)+(dukj(level))(1,rprime));
  dU+=V;

//   cerr<<"printing u_10 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (ukj(level))(q,TempukjArray); 
//     int k=1;
//     int j=0;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<< " " <<endl;
//     //  cerr<<dukj(level)(1,q)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_11 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=1;
//     int j=1;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_20 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=0;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_21 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=1;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_22 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=3;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//  q=0.5*(r+rprime);

  double UDiag=U;
  //assert(1==2);
  //  return;
  //  if (q>ukj(level).grid->End){
  //    q=ukj(level).grid->End;
  //  }
  if (s > 0.0 && q<ukj(level).grid->End) { // && q<ukj(level).grid->End){
//     if (fabs(z)>2*q){
//       z=2*q*fabs(z)/z;
//     }
//     if (s>2*q){
//       s=2*q;
//     }


    double zsquared=z*z;
    double ssquared=s*s; 
    double ssquaredinverse=1.0/ssquared;
    double Sto2k=ssquared;
    
    (ukj(level))(q,TempukjArray); 
    (dukj(level))(q,TempdukjArray); 
    ////HACK! 
    //    n=0;
    for (int k=1;k<=NumTerms;k++){  
     
      double Zto2j=1;
      double currS=Sto2k;
      
      for (int j=0;j<=k;j++){
	// indexing into the 2darray
	double Ucof  = TempukjArray(k*(k+1)/2+j+1);
	
	double dUcof = TempdukjArray(k*(k+1)/2+j+1);
	U+=(Ucof)*Zto2j*currS;
	dU+=(dUcof)*Zto2j*currS; //+V = HACK!
        //cerr << "dU " << dU << " dUcof " << dUcof << " U " << U << " Ucof " << Ucof << " Zto2j " << Zto2j << " currS " << currS << " k " << k << " j " << j << endl;
	Zto2j*=zsquared;
	currS=currS*ssquaredinverse;				
      }				
      Sto2k=Sto2k*ssquared;
    } 
  } 
//   cerr<<U<<" "<<UDiag<<" "<<U-UDiag<<endl;
//   if ((U-UDiag)>1e-1){
//     cerr<<"ERROR! ERROR!"<<endl;
//     cerr<<dU<<" "<<V<<" "<<r<<" "<<rprime<<" "<<s<<" "<<z<<" "<<q<<endl;
//   }
}

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void DavidPAClass::calcUsqzFast(double s,double q,double z,int level,
				double &U)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  U=0.0;
  double r=q+0.5*z;
  double rprime=q-0.5*z;
  if (r > ukj(level).grid->End) {
    r=ukj(level).grid->End;
  }
  if (rprime > ukj(level).grid->End) {
    rprime=ukj(level).grid->End;
  }
  // This is the endpoint action 
  if ((rprime < rmin) || (r < rmin)){
    //    cerr<<"I'm less then the min. Maybe this is messing me up\n";
    //U = 5000.0;
    U = ukj(level)(1,rmin);
    //cerr << "DavidPAClass: rprime < rmin; assigning U_0 = " << U << endl;
    return;
  }
  U+= 0.5*(ukj(level)(1,r)+ukj(level)(1,rprime));
  double UDiag=U;
  if (s > 0.0 && q<ukj(level).grid->End) { // && q<ukj(level).grid->End){
    double zsquared=z*z;
    double ssquared=s*s; 
    double ssquaredinverse=1.0/ssquared;
    double Sto2k=ssquared;
    (ukj(level))(q,TempukjArray); 
    for (int k=1;k<=NumTerms;k++){  
      double Zto2j=1;
      double currS=Sto2k;
      for (int j=0;j<=k;j++){
	// indexing into the 2darray
	double Ucof  = TempukjArray(k*(k+1)/2+j+1);
	U+=(Ucof)*Zto2j*currS;
	Zto2j*=zsquared;
	currS=currS*ssquaredinverse;				
      }				
      Sto2k=Sto2k*ssquared;
    } 
  } 
}


double 
DavidPAClass::UDiag_exact(double q,int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  double U=0.0;
  double r=q;
  if (r > ukj(level).grid->End) {
    r=ukj(level).grid->End;
  }
  // This is the endpoint action 
  if (r < rmin){
    U = 5000.0;
    return U;
  }
  U+= (ukj(level)(1,r));
  return U;
}



double
DavidPAClass::Udiag (double q, int level)
{
  level=level+TauPos;
  // This is the endpoint action   
  if (q < UdiagSpline(level).Start())
    return UdiagSpline(level)(UdiagSpline(level).Start());
  else if (q > UdiagSpline(level).End())
    return UdiagSpline(level)(UdiagSpline(level).End());
  else
    return UdiagSpline(level)(q);
}

double
DavidPAClass::dUdiag_fast (double q, int level)
{
  level=level+TauPos;
  // This is the endpoint action   
  if (q < dUdiagSpline(level).Start()) 
    return dUdiagSpline(level)(dUdiagSpline(level).Start());
  else if (q > dUdiagSpline(level).End())
    return dUdiagSpline(level)(dUdiagSpline(level).End());
  else
    return dUdiagSpline(level)(q);
}



void DavidPAClass::ReadSamplingTable(string fileName)
{
  verr<<"Reading the sampling table"<<endl;
  IOSectionClass in;
  bool success=in.OpenFile(fileName.c_str());
  if (!success){
    cerr<<"Expected but can not find a pairaction file";
    assert(1==2);
  }
  int startLevel=-1;
  assert(in.ReadVar("SamplingTau",SamplingTau));
  for (int i=0;i<SamplingTau.size();i++){
    //cerr<<"Tau diffs: "<<SamplingTau(i)<<" "<<DesiredTau<<" "<<SamplingTau(i)-DesiredTau<<endl;
    if (fabs(SamplingTau(i)-DesiredTau)<1e-4){
      verr<<"The sampling tau I've chosen is "<<SamplingTau(i);
      startLevel=i;
    }
  }
  assert(startLevel!=-1);
  in.OpenSection("dUdR");
  //Now we want to read the grid
  assert(in.OpenSection("Grid"));
  string gridType;
  double gridStart;
  double gridEnd;
  int numGridPoints;
  assert(in.ReadVar("Type",gridType));
  assert(gridType=="Linear");
  assert(in.ReadVar("start",gridStart));
  assert(in.ReadVar("end",gridEnd));
  assert(in.ReadVar("NumPoints",numGridPoints));
  LinearGrid *theGrid;
  theGrid=new LinearGrid(gridStart,gridEnd,numGridPoints);
  in.CloseSection();

  Array<double,2> tempdudR;
  Array<double,2> tempdudR_movers;
  in.ReadVar("dUdR",tempdudR);
  in.ReadVar("dUdR_movers",tempdudR_movers);
  //We should now assert that the tau array has the correct number
  //of tau's in it.
  assert(tempdudR.extent(0)==NumTau);
  assert(tempdudR_movers.extent(0)==NumTau);
  dUdRTimesSigmaSpline.resize(NumTau-startLevel);
  dUdRTimesSigmaSpline_movers.resize(NumTau-startLevel);
  for (int level=startLevel;level<NumTau;level++){
    dUdRTimesSigmaSpline(level-startLevel).Init(theGrid,tempdudR(level,Range::all()));
    dUdRTimesSigmaSpline_movers(level-startLevel).Init(theGrid,tempdudR_movers(level,Range::all()));
  }
  ///We've read the dUdR and now will read the d2UdR2



  Array<double,2> tempd2UdR2;
  Array<double,2> tempd2UdR2_movers;
  in.ReadVar("d2UdR2_nonmovers",tempd2UdR2);
  in.ReadVar("d2UdR2_movers",tempd2UdR2_movers);
  //We should now assert that the tau array has the correct number
  //of tau's in it.
  assert(tempd2UdR2.extent(0)==NumTau);
  assert(tempd2UdR2_movers.extent(0)==NumTau);
  d2UdR2TimesSigmaSpline.resize(NumTau-startLevel);
  d2UdR2TimesSigmaSpline_movers.resize(NumTau-startLevel);
  for (int level=startLevel;level<NumTau;level++){
    d2UdR2TimesSigmaSpline(level-startLevel).Init(theGrid,tempd2UdR2(level,Range::all()));
    d2UdR2TimesSigmaSpline_movers(level-startLevel).Init(theGrid,tempd2UdR2_movers(level,Range::all()));
  }
  
  


  
  in.CloseSection();
  in.CloseFile();
  verr<<"Left the sampling table"<<endl;
}

void DavidPAClass::ReadDavidSquarerFile(string DMFile)
{
  ///  cerr<<"I AM READING DAVID SQUARER FILE"<<endl;
  double tau; //used to be in the base clase
  double smallestTau;
  ifstream infile;
  //cout <<DMFile<<endl;
  infile.open(DMFile.c_str());  
  if (infile.fail()){
    cerr<<"CAN'T OPEN THE FILE!!!!!!!!!!";
  }
  
  string numOfFitsString=SkipTo(infile,"SQUARER");
  //  cerr<<GetNextWord(numOfFitsString)<<endl;
  GetNextWord(numOfFitsString);
  //  double LowTau=atof(GetNextWord(numOfFitsString));
  ///  cerr<<GetNextWord(numOfFitsString)<<endl;
  GetNextWord(numOfFitsString);
  //  cerr<<LowTau<<endl;
  //  cerr<<GetNextWord(numOfFitsString)<<endl;
  GetNextWord(numOfFitsString);
  //  cerr<<GetNextWord(numOfFitsString)<<endl;
  GetNextWord(numOfFitsString);

  int numOfFits=GetNextInt(numOfFitsString);
  n = numOfFits;
  // Read in  the potential
  Array<double,1> potential;
  string potGridString = SkipTo(infile, "RANK");
  GetNextWord(potGridString);
  GetNextWord(potGridString);//HACK?
  int numPotPoints = GetNextInt(potGridString);
  potential.resize(numPotPoints);

  SkipTo(infile, "potential");

  for (int i=0; i<numPotPoints; i++){
    infile >> potential(i);

  }

  //  string NDERIVString = SkipTo(infile,"NDERIV");


  //  NDERIVString.erase(NDERIVString.find("NDERIV"),strlen("NDERIV"));

  ///  2*(NDERIV+1);
  Grid *theGrid;

  for (int counter=0;counter<=numOfFits;counter++){ //Get the U's 
    string RankString =SkipTo(infile,"RANK");
    int theRank=GetNextInt(RankString);
    //cout<<theRank<<endl;

    if (theRank!=3){
      //cerr<<"ERROR! ERROR! Rank was not 3" << endl;
      counter--;
    }
    else {
      int NumGridPoints=GetNextInt(RankString);
      int NumUKJ=GetNextInt(RankString);
      NumTau=GetNextInt(RankString);

      
      
      string RGridString =SkipTo(infile,"GRID 1");
      string GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      double startGrid = GetNextDouble(RGridString);
      double endGrid = GetNextDouble(RGridString);
      GridType.resize(3);
      if (GridType=="LIN"){
	theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
      }
      else if (GridType == "LOG"){
	//cout<<"We're really in log grid here\n";
	double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
	//cerr << "delta = " << delta << endl;
	theGrid = new LogGrid(startGrid,delta,NumGridPoints);
      }
      else {
	cerr << "Unrecognized grid type in ReadDavidSquarerFile (text).\n";
	cerr << "GridType = \"" << GridType << "\"\n";
        abort();
      }
	  
      
      string TauGridString = SkipTo(infile,"GRID   3"); //We hope this is a log grid
      GetNextWord(TauGridString);
      GetNextWord(TauGridString); /// takes out the Grid  3
      string shouldBeLog;
      if  ((shouldBeLog=GetNextWord(TauGridString))!="LOG"){
	cerr<<"ERROR!!! ERROR!!! The tau grid is not a LOG Grid\n";
	cerr<<shouldBeLog<<endl;
      }
      smallestTau=GetNextDouble(TauGridString);
      double largestTau=GetNextDouble(TauGridString);
      int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
      if (NumTau!=numTauCalc){
	
	cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
	cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
      }
      string beginString=SkipTo(infile,"BEGIN density matrix table");
      int NMax=GetNextInt(beginString); //This is magically the most accurate fit i.e. NDERIV-1
      if (GetNextInt(beginString)!=1){ //i.e. if it's not U
	cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";
      }
      Array<double,3> tempUkj(NumGridPoints,NumUKJ,NumTau);
      //      cerr<<"NumTau is"<<NumTau<<endl;
      ukj.resize(NumTau);
      UdiagSpline.resize(NumTau);
      dUdiagSpline.resize(NumTau);
      ////???      ukj.resize(NumUKJ+1);
      ReadFORTRAN3Tensor(infile,tempUkj);
      Array<double,3> tempUkj2(NumGridPoints,NumUKJ+1,NumTau);
      for(int i=0; i<NumTau; i++){
	tempUkj2(Range::all(),0,i) = potential;
      }
      tempUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempUkj;
      tempUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0;
      Array<double,1> startDeriv(NumUKJ+1); //I think this is the right number of grid poins
      startDeriv=5.0e30;
      Array<double,1> endDeriv(NumUKJ+1);
      endDeriv=0.0;

      ///////      tau=largestTau; //HACK!
      tau=smallestTau;
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){//the -3 here is a HACK!
	if (NMax==2){ //MORE HACK!
	  //ukj(levelCounter).-Init(theGrid,tempUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
	  ukj(levelCounter).Init(theGrid,tempUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
	}
	tau=tau*2; //HACK!
      }
      //      tau=smallestTau; HACK REMOVAL!
      n=NMax;
      
    }
  }



  for (int counter=0;counter<=numOfFits;counter++){ //Get the beta derivative of U's 
    
    string RankString =SkipTo(infile,"RANK");
    int theRank=GetNextInt(RankString);
    //cout<<theRank<<endl;

    if (theRank!=3){
      //cerr<<"ERROR! ERROR! Rank was not 3" << endl;
      counter--;
    }
    else {
      int NumGridPoints=GetNextInt(RankString);
      int NumUKJ=GetNextInt(RankString);
      NumTau=GetNextInt(RankString);
      
      
      string RGridString =SkipTo(infile,"GRID 1");
      string GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      double startGrid = GetNextDouble(RGridString);
      double endGrid = GetNextDouble(RGridString);
      GridType.resize(3);
      if (GridType=="LIN"){
	theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
      }
      else if (GridType == "LOG"){
	//cout<<"We're really in log grid here\n";
	double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
	//cerr << "delta = " << delta << endl;
	theGrid = new LogGrid(startGrid,delta,NumGridPoints);
      }
      else {
	cerr << "Unrecognized grid type in ReadDavidSquarerFile.\n";
	cerr << "GridType = \"" << GridType << "\"\n";
        abort();
      }
	  
      
      string TauGridString = SkipTo(infile,"GRID   3"); //We hope this is a log grid
      GetNextWord(TauGridString);
      GetNextWord(TauGridString); /// takes out the Grid  3
      string shouldBeLog;
      if  ((shouldBeLog=GetNextWord(TauGridString))!="LOG"){
	cerr<<"ERROR!!! ERROR!!! The tau grid is not a LOG Grid\n";
	cerr<<shouldBeLog<<endl;
      }
      smallestTau=GetNextDouble(TauGridString);
      double largestTau=GetNextDouble(TauGridString);
      int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
      ////      cerr<<"The largest and smallest tau are respectively"
      ////	  <<largestTau<<" "<<smallestTau<<endl;
      if (NumTau!=numTauCalc){
	
	cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
	cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
      }
      string beginString=SkipTo(infile,"BEGIN density matrix table");
      int NMax=GetNextInt(beginString); //This is magically the most accurate fit i.e. NDERIV-1
      ///      cerr<<"NMAX is "<<NMax<<endl;
      if (GetNextInt(beginString)!=2){ //i.e. if it's not U
	cerr<<"ERROR!!! ERROR!!! We didn't get the beta derivative.\n";
      }
      Array<double,3> tempdUkj(NumGridPoints,NumUKJ,NumTau);
      Array<double,3> tempdUkj2(NumGridPoints,NumUKJ+1,NumTau);
      if (NMax==2){
	TempukjArray.resize(NumUKJ+1);      
	TempdukjArray.resize(NumUKJ+1);      
      }
      dukj.resize(NumTau);
      ///???dukj.resize(NumUKJ+1);
      Array<double,1> startDeriv(NumUKJ+1); //I think this is the right number of grid poins
      startDeriv=5.0e30;
      Array<double,1> endDeriv(NumUKJ+1);
      endDeriv=0.0;

      ReadFORTRAN3Tensor(infile,tempdUkj);
      /////      tau=largestTau; //HACK
      tau=smallestTau;
      for(int i=0; i<NumTau; i++){ //HACK!
	tempdUkj2(Range::all(),0,i) = potential;
	///	cerr<<"Current tau is "<<tau<<" "<<i<<endl;
	if (fabs(tau-DesiredTau)<1e-6){
	  ///	  cerr<<"The tau I've chosen is "<<tau;
	  TauPos=i;
	}
	tau=tau*2; //HACK!
      }
      ///      cerr<<"I'm about ot actually initialize dukj now!"<<endl;
      tempdUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempdUkj;
      tempdUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0; ///NOT SURE ABOUT THIS!!!
      const int numDiagPoints = 2000;
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
	
	if (NMax==2){ //NMax again
	  dukj(levelCounter).Init(theGrid,tempdUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
	}
      }
      //      tau=smallestTau; HACK REMOVAL!
      ////      cerr<<"My tau is "<<tau<<endl;
      n=NMax;
      
    }
    
  }
  Potential.resize(potential.size());
  for (int counter=0;counter<potential.size();counter++){
    Potential(counter)=potential(counter);
  }
  tau=smallestTau;
  for (int i=0;i<TauPos;i++){
    tau *= 2;
  }
  ///  cerr << "NumTau = " << NumTau << endl;
  for (int level=0; level<NumTau; level++) {
    const int numDiagPoints = 20000;
    Array<double,1> udiag(numDiagPoints);
    Array<double,1> dUdiag(numDiagPoints);
    double start = ukj(level).grid->Start;
    double end   = ukj(level).grid->End;
    double dr = (end-start)/(double)(numDiagPoints-1);
    for (int j=0; j<numDiagPoints; j++) {
      double r = start + (double)j * dr;
      calcUsqzFast (0.0, r, 0.0, level-TauPos, udiag(j));
      dUdiag(j)=dU(r,0.0,0.0,level-TauPos);
    }
    UdiagSpline(level).Init (start, end, udiag);
    dUdiagSpline(level).Init (start, end, dUdiag);
  }
  verr<<"I've selected a tau of "<<tau<< "in the PairAction file"<<endl;
  ///  cerr<<"TauPos is "<<TauPos<<endl;
}

void DavidPAClass::ReadLongRangeHDF5(IOSectionClass &in)
{



  assert(in.ReadVar("kPoints",kVals));  
  kCutoff=0.0;
  if (!in.ReadVar("kcut",kCutoff))
    cerr<<"kCutoff is not specified.  Files could be inconsistent";
  assert(in.ReadVar("u_k",uk_long));
  assert(in.ReadVar("Box",LongRangeBox));
  assert(in.ReadVar("mass1",LongRangeMass1));
  assert(in.ReadVar("mass2",LongRangeMass2));
  assert(in.ReadVar("ndim",LongRangeDim));
}



void DavidPAClass::ReadDavidSquarerFileHDF5(string DMFile)
{
  TauPos = -1;
  //cerr << "READDAVIDSQUARERFILE -- HDF5 -- " << endl;
  double tau; //used to be in the base clase
  double smallestTau;
  IOSectionClass in;
  if(!in.OpenFile(DMFile.c_str())) {
    cerr << "ERROR: Could not find pair action file " << DMFile << ". Aborting..." << endl;
    abort();
  }
  int numOfFits;
  assert(in.OpenSection("Squarer"));
  assert(in.ReadVar("NumFits", numOfFits));
  Vimage=0.0;
  if (!in.ReadVar("vimage",Vimage))
    cerr<<"WARNING! This is a pre-vimage dm matrix and will not give the right correction to the energy"<<endl;
  in.CloseSection();
  n = numOfFits;

  // Read in  the potential
  Array<double,1> potential;
  assert(in.OpenSection("Potential"));
  assert(in.ReadVar("Data", potential));
  in.CloseSection();

  // Get the U's
  Grid *theGrid;
  Array<double,1> Taus;
  for (int counter=0;counter<=numOfFits;counter++) {
    ostringstream stream;
    stream << "Ukj" << counter;
    string SectionTitle(stream.str());
    assert(in.OpenSection(SectionTitle));

    int theRank;
    in.ReadVar("Rank",theRank);
    if (theRank!=3){
      cerr<<"ERROR! ERROR! Rank was not 3" << endl;
      counter--;
    }
    else {
      int NumGridPoints, NumUKJ;
      in.ReadVar("NumUkj",NumUKJ);
      in.ReadVar("NumTau",NumTau);

      assert(in.OpenSection("Grid"));
      string GridType;
      double startGrid, endGrid;
      in.ReadVar("NumGridPoints",NumGridPoints);
      in.ReadVar("Type",GridType);
      in.ReadVar("Start",startGrid);
      in.ReadVar("End",endGrid);
      in.CloseSection();
      GridType.resize(3);
      if (GridType=="LIN") {
        theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
      }
      else if (GridType == "LOG") {
        double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
        theGrid = new LogGrid(startGrid,delta,NumGridPoints);
      }
      else {
        cerr << "Unrecognized grid type in ReadDavidSquarerFile. (u hdf5)\n";
        cerr << "GridType = \"" << GridType << "\"\n";
        abort();
      }

      in.ReadVar("Taus",Taus);
      smallestTau=Taus(0);
      double largestTau=Taus(Taus.size()-1);
      int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
      if (NumTau!=numTauCalc) {
        cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
        cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
      }
      int NMax;
      in.ReadVar("NMax",NMax);
      int derv;
      in.ReadVar("Derv",derv);
      if (derv!=1){ //i.e. if it's not U
        cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";
      }
      Array<double,3> tempUkj(NumGridPoints,NumUKJ,NumTau);
      ukj.resize(NumTau);
      UdiagSpline.resize(NumTau);
      dUdiagSpline.resize(NumTau);
      in.ReadVar("Data",tempUkj);
      Array<double,3> tempUkj2(NumGridPoints,NumUKJ+1,NumTau);
      for(int i=0; i<NumTau; i++) {
        tempUkj2(Range::all(),0,i) = potential;
      }
      tempUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempUkj;
      tempUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0;
      Array<double,1> startDeriv(NumUKJ+1); //I think this is the right number of grid points
      startDeriv=5.0e30;
      Array<double,1> endDeriv(NumUKJ+1);
      endDeriv=0.0;

      tau=smallestTau;
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){//the -3 here is a HACK!
        if (NMax==2) { //MORE HACK!
          ukj(levelCounter).Init(theGrid,tempUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
        }
        tau=tau*2; //HACK!
      }
      //      tau=smallestTau; HACK REMOVAL!
      n=NMax;

    }
    in.CloseSection();
  }

  //Get the beta derivative of U's
  //cerr << "Reading du" << endl;
  //int checkNumdU = in.CountSections("dUkj_dBeta");
  //assert(checkNumdU == (numOfFits+1));
  for (int counter=0;counter<=numOfFits;counter++){
    //cerr << "  " << counter << " of " << numOfFits << endl;
    ostringstream stream;
    stream << "dUkjdBeta" << counter;
    string SectionTitle(stream.str());
    //cerr << "Opening section " << SectionTitle << "||" << endl;
    assert(in.OpenSection(SectionTitle));

    int theRank;
    //cerr<<"Reading du rank"<<endl;
    in.ReadVar("Rank",theRank);
    //cerr<<"du rank read"<<endl;
    assert(theRank == 3);
    int NumGridPoints, NumUKJ;
    in.ReadVar("NumUkj",NumUKJ);
    in.ReadVar("NumTau",NumTau);

    assert(in.OpenSection("Grid"));
    string GridType;
    double startGrid, endGrid;
    in.ReadVar("NumGridPoints",NumGridPoints);
    in.ReadVar("Type",GridType);
    in.ReadVar("Start",startGrid);
    in.ReadVar("End",endGrid);
    in.CloseSection();
    GridType.resize(3);
    if (GridType=="LIN"){
      theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
    }
    else if (GridType == "LOG"){
      double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
      theGrid = new LogGrid(startGrid,delta,NumGridPoints);
    }
    else {
      cerr << "Unrecognized grid type in ReadDavidSquarerFile (du hdf5).\n";
      cerr << "GridType = \"" << GridType << "\"\n";
      abort();
    }
    //cerr<<"Reading taus"<<endl;
    in.ReadVar("Taus",Taus);
    //cerr<<"taus rad"<<endl;
    smallestTau=Taus(0);
    double largestTau=Taus(Taus.size()-1);
    int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
    if (NumTau!=numTauCalc){
      cerr <<"ERROR!!! ERROR!!! num tau inconsistency \n";
      cerr <<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
    }
    int NMax;
    in.ReadVar("NMax",NMax);
    int derv;
    in.ReadVar("Derv",derv);
    if (derv!=2){ //i.e. if it's not U
      cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";
    }
    Array<double,3> tempdUkj(NumGridPoints,NumUKJ,NumTau);
    Array<double,3> tempdUkj2(NumGridPoints,NumUKJ+1,NumTau);
    if (NMax==2){
      TempukjArray.resize(NumUKJ+1);
      TempdukjArray.resize(NumUKJ+1);
    }
    dukj.resize(NumTau);
    Array<double,1> startDeriv(NumUKJ+1); //I think this is the right number of grid poins
    startDeriv=5.0e30;
    Array<double,1> endDeriv(NumUKJ+1);
    endDeriv=0.0;
    in.ReadVar("Data",tempdUkj);

    tau=smallestTau;
    for(int i=0; i<NumTau; i++){ //HACK!
      tempdUkj2(Range::all(),0,i) = potential;
      if (fabs(tau-DesiredTau)<1e-6){
        TauPos=i;
      }
      tau=tau*2; //HACK!
    }
    tempdUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempdUkj;
    tempdUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0; ///NOT SURE ABOUT THIS!!!
    const int numDiagPoints = 20000;
    for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
      if(NMax==2){ //NMax again
        dukj(levelCounter).Init(theGrid,tempdUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
      }
    }
    n=NMax;
    in.CloseSection();
  }

  Potential.resize(potential.size());
  for (int counter=0;counter<potential.size();counter++){
    Potential(counter)=potential(counter);
  }
  tau=smallestTau;
  if(TauPos == -1) {
    cerr << "Tau of " << DesiredTau << " not found.  Possibilities are " << Taus << endl;
    cerr << "ABORTING" << endl;
    exit(1);
  }
  for (int i=0;i<TauPos;i++){
    tau *= 2;
  }
  //cerr<<"pre splining"<<endl;
  //cerr << "NumTau = " << NumTau << endl;
  for (int level=0; level<NumTau; level++) {
    //cerr<<"level is "<<level<<endl;
    const int numDiagPoints = 20000;
    Array<double,1> udiag(numDiagPoints);
    Array<double,1> dUdiag(numDiagPoints);
    double start = ukj(level).grid->Start;
    double end   = ukj(level).grid->End;
    double dr = (end-start)/(double)(numDiagPoints-1);
    for (int j=0; j<numDiagPoints; j++) {
      double r = start + (double)j * dr;
      calcUsqzFast (0.0, r, 0.0, level-TauPos, udiag(j));
      dUdiag(j)=dU(r,0.0,0.0,level-TauPos);
    }
    UdiagSpline(level).Init (start, end, udiag);
    dUdiagSpline(level).Init (start, end, dUdiag);
    //cerr<<"Bottom of loop"<<endl;
  }
  //cerr<<"done with loop"<<endl;
  //verr<<"I've selected a tau of "<<tau<< "in the PairAction file"<<endl;
  //cerr<<"TauPos is "<<TauPos<<endl;
  if (in.OpenSection("LongRange")){
    HasLongRange=true;
    //cerr<<"In long range finding"<<endl;
    ReadLongRangeHDF5(in);
    in.CloseSection();
    //cerr<<"done with long range finding"<<endl;
  }
  else
    HasLongRange=false;
  //cerr << "Leaving DavidSquarer HDF5 read" << endl;
  //

  /// Print out some action values
  // PrintVals(0.0,3.0,0.1)
}

void DavidPAClass::PrintVals(double begin, double end, double dx)
{
  double x = begin;
  while (x<=end) {
    cout << x << " " << Udiag(x,0) << " " << dUdiag(x,0) << " " << dU(x,0.0,0.0,0) << endl;
    x += dx;
  }
}

// double DavidPAClass::Udiag(double q, int level)
// {
//   level=level+TauPos;
//   double rmin = ukj(level).grid->Start;
  
//   if (q > ukj(level).grid->End)
//     q = ukj(level).grid->End;

//   // This is the endpoint action   
//   if (q < rmin) 
//     return 5000.0;
  
//   return ukj(level)(1,q); 
// }

double DavidPAClass::Udiag_p(double q, int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  
  if (q > ukj(level).grid->End)
    q = ukj(level).grid->End;

  // This is the endpoint action   
  if (q < rmin) 
    //return 5000.0;
    q = rmin;
  
  return ukj(level).Deriv(1,q); 
}

double DavidPAClass::Udiag_pp(double q, int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  
  if (q > ukj(level).grid->End)
    q = ukj(level).grid->End;

  // This is the endpoint action   
  if (q < rmin) 
    //return 5000.0;
    q = rmin;
  
  return ukj(level).Deriv2(1,q); 
}


double DavidPAClass::dUdiag(double q, int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  
  if (q > ukj(level).grid->End)
    q = ukj(level).grid->End;

  if (q < rmin) 
    //return 5000.0;
    q = rmin;
  
  // This is the endpoint action
  //                            this is V(q)
  //                             |
  //                             v   
  return dukj(level)(1,q) + ukj(level)(0,q);
}


double DavidPAClass::dUdiag_p(double q, int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  
  if (q > ukj(level).grid->End)
    q = ukj(level).grid->End;

  if (q < rmin) 
    //return 5000.0;
    q = rmin;
  
  // This is the endpoint action
  //                            this is V(q)
  //                             |
  //                             v   
  return dukj(level).Deriv(1,q) + ukj(level).Deriv(0,q);
}


double DavidPAClass::dUdiag_pp(double q, int level)
{
  level=level+TauPos;
  double rmin = ukj(level).grid->Start;
  
  if (q > ukj(level).grid->End)
    q = ukj(level).grid->End;

  if (q < rmin) 
    //return 5000.0;
    q = rmin;
  
  // This is the endpoint action
  //                            this is V(q)
  //                             |
  //                             v   
  return dukj(level).Deriv2(1,q) + ukj(level).Deriv2(0,q);
}

