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
  calcUsqzFast(s,q,z,level,uTemp);
  return uTemp;
}


double DavidPAClass::V(double r)
{
  double rMin = ukj(0).grid->Start;
  double rMax = ukj(0).grid->End;
  if (r>rMax)
    r = rMax;
  if (r<rMin)
    r = rMin;
  return ukj(0)(0,r);
}


double DavidPAClass::dU(double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp, duTemp, vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  return duTemp;
}


bool DavidPAClass::IsLongRange()
{
  return false;
}


void DavidPAClass::Derivs (double q, double z, double s2, int level, double &d_dq, double &d_dz)
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
    for (int k=1; k<=nTerms; k++) {  
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


void DavidPAClass::DerivsFD (double q, double z, double s2, int level, double &d_dq, double &d_dz)
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


void DavidPAClass::Derivs (double q, double z, double s2, int level, double &d_dq, double &d_dz, double &d_ds)
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
    for (int k=1; k<=nTerms; k++) {  
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


void DavidPAClass::DerivsFD (double q, double z, double s2, int level, double &d_dq, double &d_dz, double &d_ds)
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
void DavidPAClass::calcUsqz(double s, double q, double z, int level, double &U, double &dU, double &V)
{
  level = level+TauPos;
  double r = q + 0.5*z;
  double rprime = q - 0.5*z;

  // Limits
  double rMax = ukj(level).grid->End;
  if (r > rMax)
    r = rMax;
  if (rprime > rMax)
    rprime = rMax;
  double rmin = ukj(level).grid->Start;
  if (rprime < rmin)
    rprime = rmin;
  if(r < rmin)
    r = rmin;

  // This is the endpoint action
  V = 0.5*(ukj(level)(0,r) + ukj(level)(0,rprime));
  U = 0.5*(ukj(level)(1,r) + ukj(level)(1,rprime));
  dU = 0.5*((dukj(level))(1,r) + (dukj(level))(1,rprime));

  // Compensate for potential, which is subtracted from diaganal action in dm file.
  dU += V;

  // Add in off-diagonal terms
  if (s > 0.0 && q<rMax) {
    double zsquared = z*z;
    double ssquared = s*s;
    double ssquaredinverse = 1./ssquared;
    double Sto2k = ssquared;
    (ukj(level))(q,TempukjArray);
    (dukj(level))(q,TempdukjArray);
    for (int k=1; k<=nTerms; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (int j=0; j<=k; j++){
        // indexing into the 2darray
        double Ucof = TempukjArray(k*(k+1)/2+j+1);
        double dUcof = TempdukjArray(k*(k+1)/2+j+1);
        U += (Ucof)*Zto2j*currS;
        dU += (dUcof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}


/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void DavidPAClass::calcUsqzFast(double s, double q, double z, int level, double &U)
{
  level = level+TauPos;
  double r = q + 0.5*z;
  double rprime = q - 0.5*z;

  // Limits
  if (r > ukj(level).grid->End)
    r = ukj(level).grid->End;
  double rMax = ukj(level).grid->End;
  if (rprime > rMax)
    rprime = rMax;
  double rMin = ukj(level).grid->Start;
  if (rprime < rMin)
    rprime = rMin;
  if(r < rMin)
    r = rMin;

  // This is the endpoint action
  U = 0.5*(ukj(level)(1,r) + ukj(level)(1,rprime));

  // Add in off-diagonal terms
  if (s>0.0 && q<rMax) { // && q<ukj(level).grid->End){
    double zsquared = z*z;
    double ssquared = s*s;
    double ssquaredinverse = 1./ssquared;
    double Sto2k = ssquared;
    (ukj(level))(q,TempukjArray);
    for (int k=1; k<=nTerms; k++) {
      double Zto2j = 1;
      double currS = Sto2k;
      for (int j=0; j<=k; j++) {
        // indexing into the 2darray
        double Ucof = TempukjArray(k*(k+1)/2 + (j+1));
        U += (Ucof)*Zto2j*currS;
        Zto2j *= zsquared;
        currS *= ssquaredinverse;
      }
      Sto2k *= ssquared;
    }
  }

}


double DavidPAClass::UDiag_exact(double q,int level)
{
  level = level+TauPos;
  double rMin = ukj(level).grid->Start;
  double rMax = ukj(level).grid->End;
  double r = q;
  if (r > rMax)
    r = rMax;
  else if (r < rMin)
    r = rMin;
  double U = (ukj(level)(1,r));
  return U;
}


double DavidPAClass::Udiag (double q, int level)
{
  level = level+TauPos;
  if (q < UdiagSpline(level).Start())
    return UdiagSpline(level)(UdiagSpline(level).Start());
  else if (q > UdiagSpline(level).End())
    return UdiagSpline(level)(UdiagSpline(level).End());
  else
    return UdiagSpline(level)(q);
}


double DavidPAClass::dUdiag_fast (double q, int level)
{
  level=level+TauPos;
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
  //We should now assert that the tau array has the correct number of tau's in it.
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
  
  string nFitsString=SkipTo(infile,"SQUARER");
  //  cerr<<GetNextWord(nFitsString)<<endl;
  GetNextWord(nFitsString);
  //  double LowTau=atof(GetNextWord(nFitsString));
  ///  cerr<<GetNextWord(nFitsString)<<endl;
  GetNextWord(nFitsString);
  //  cerr<<LowTau<<endl;
  //  cerr<<GetNextWord(nFitsString)<<endl;
  GetNextWord(nFitsString);
  //  cerr<<GetNextWord(nFitsString)<<endl;
  GetNextWord(nFitsString);

  int nFits=GetNextInt(nFitsString);
  n = nFits;
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

  for (int counter=0;counter<=nFits;counter++){ //Get the U's 
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



  for (int counter=0;counter<=nFits;counter++){ //Get the beta derivative of U's 
    
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
  kCutoff = 0.0;
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

  // Read in Squarer parameters
  IOSectionClass in;
  if(!in.OpenFile(DMFile.c_str())) {
    cerr << "ERROR: Could not find pair action file " << DMFile << ". Aborting..." << endl;
    abort();
  }
  int nFits;
  assert(in.OpenSection("Squarer"));
  assert(in.ReadVar("NumFits", nFits));
  if (nFits < nTerms) {
    cerr << "ERROR: Asking for a higher pair action expansion than available, nFits: " << nFits << " nTerms: " << nTerms << endl;
    abort();
  }
  Vimage = 0.0;
  if (!in.ReadVar("vimage",Vimage))
    cerr<<"WARNING! This is a pre-vimage dm matrix and will not give the right correction to the energy"<<endl;
  in.CloseSection();

  // Read in the potential
  Array<double,1> potential;
  assert(in.OpenSection("Potential"));
  assert(in.ReadVar("Data", potential));
  in.CloseSection();

  // Get the U's
  ostringstream stream;
  stream << "Ukj" << nTerms;
  string SectionTitle(stream.str());
  assert(in.OpenSection(SectionTitle));

  // U Parameters
  int theRank;
  in.ReadVar("Rank",theRank);
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
  Grid *UGrid;
  if (GridType=="LIN") {
    UGrid = new LinearGrid(startGrid,endGrid,NumGridPoints);
  }
  else if (GridType == "LOG") {
    double delta = pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
    UGrid = new LogGrid(startGrid,delta,NumGridPoints);
  }
  else {
    cerr << "Unrecognized grid type in ReadDavidSquarerFile. (u hdf5)\n";
    cerr << "GridType = \"" << GridType << "\"\n";
    abort();
  }

  // Make sure tau is correct
  Array<double,1> Taus;
  in.ReadVar("Taus",Taus);
  double smallestTau = Taus(0);
  double largestTau = Taus(Taus.size()-1);
  int numTauCalc = (int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0);
  if (NumTau!=numTauCalc) {
    cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
    cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
  }
  TauPos = -1;
  double tau = smallestTau;
  for(int i=0; i<NumTau; i++) {
    if (fabs(tau-DesiredTau)<1e-6)
      TauPos = i;
    tau *= 2; //HACK!
  }
  if(TauPos == -1) {
    cerr << "Tau of " << DesiredTau << " not found.  Possibilities are " << Taus << endl;
    abort();
  }

  int NMax;
  in.ReadVar("NMax",NMax);
  int derv;
  in.ReadVar("Derv",derv);
  if (derv!=1) //i.e. if it's not U
    cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";

  // Initiate U Splines
  Array<double,3> tempUkj(NumGridPoints,NumUKJ,NumTau);
  ukj.resize(NumTau);
  UdiagSpline.resize(NumTau);
  dUdiagSpline.resize(NumTau);
  in.ReadVar("Data",tempUkj);
  Array<double,3> tempUkj2(NumGridPoints,NumUKJ+1,NumTau);
  for(int i=0; i<NumTau; i++)
    tempUkj2(Range::all(),0,i) = potential;
  tempUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempUkj;
  tempUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0;
  Array<double,1> startDeriv(NumUKJ+1); //I think this is the right number of grid points
  startDeriv=5.0e30;
  Array<double,1> endDeriv(NumUKJ+1);
  endDeriv=0.0;
  for (int levelCounter=0; levelCounter<NumTau; levelCounter++)
    ukj(levelCounter).Init(UGrid,tempUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
  in.CloseSection();

  // Get the beta derivative of U's
  ostringstream dUStream;
  dUStream << "dUkjdBeta" << nTerms;
  SectionTitle = dUStream.str();
  assert(in.OpenSection(SectionTitle));

  // dU Parameters
  in.ReadVar("Rank",theRank);
  assert(theRank == 3);
  in.ReadVar("NumUkj",NumUKJ);
  in.ReadVar("NumTau",NumTau);
  assert(in.OpenSection("Grid"));
  in.ReadVar("NumGridPoints",NumGridPoints);
  in.ReadVar("Type",GridType);
  in.ReadVar("Start",startGrid);
  in.ReadVar("End",endGrid);
  in.CloseSection();
  GridType.resize(3);
  Grid *dUGrid;
  if (GridType=="LIN") {
    dUGrid = new LinearGrid(startGrid,endGrid,NumGridPoints);
  }
  else if (GridType == "LOG") {
    double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
    dUGrid = new LogGrid(startGrid,delta,NumGridPoints);
  }
  else {
    cerr << "Unrecognized grid type in ReadDavidSquarerFile (du hdf5).\n";
    cerr << "GridType = \"" << GridType << "\"\n";
    abort();
  }
  in.ReadVar("Taus",Taus);
  smallestTau = Taus(0);
  largestTau = Taus(Taus.size()-1);
  numTauCalc = (int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
  if (NumTau!=numTauCalc) {
    cerr <<"ERROR!!! ERROR!!! num tau inconsistency \n";
    cerr <<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
  }
  in.ReadVar("NMax",NMax);
  in.ReadVar("Derv",derv);
  if (derv!=2) {
    cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";
    abort();
  }
  Array<double,3> tempdUkj(NumGridPoints,NumUKJ,NumTau);
  Array<double,3> tempdUkj2(NumGridPoints,NumUKJ+1,NumTau);
  TempukjArray.resize(NumUKJ+1);
  TempdukjArray.resize(NumUKJ+1);
  dukj.resize(NumTau);
  startDeriv.resize(NumUKJ+1); //I think this is the right number of grid poins
  startDeriv=5.0e30;
  endDeriv.resize(NumUKJ+1);
  endDeriv=0.0;
  in.ReadVar("Data",tempdUkj);

  // Initiate dU Splines
  for(int i=0; i<NumTau; i++)
    tempdUkj2(Range::all(),0,i) = potential;
  tempdUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempdUkj;
  tempdUkj2(NumGridPoints-1,Range::all(),Range::all()) = 0.0; ///NOT SURE ABOUT THIS!!!
  const int numDiagPoints = 20000;
  for (int levelCounter=0;levelCounter<NumTau;levelCounter++)
    dukj(levelCounter).Init(dUGrid,tempdUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
  in.CloseSection();

  // Initiate Potential
  Potential.resize(potential.size());
  for (int counter=0; counter<potential.size(); counter++)
    Potential(counter) = potential(counter);


  // Build diagonal splines
  for (int level=0; level<NumTau; level++) {
    const int numDiagPoints = 20000;
    Array<double,1> udiag(numDiagPoints);
    Array<double,1> dUdiag(numDiagPoints);
    double start = ukj(level).grid->Start;
    double end = ukj(level).grid->End;
    double dr = (end-start)/(double)(numDiagPoints-1);
    for (int j=0; j<numDiagPoints; j++) {
      double r = start + (double)j * dr;
      calcUsqzFast(0.0,r,0.0,level-TauPos,udiag(j));
      dUdiag(j) = dU(r,0.0,0.0,level-TauPos);
    }
    UdiagSpline(level).Init(start, end, udiag);
    dUdiagSpline(level).Init(start, end, dUdiag);
  }

  // Extra Long Range stuff
  if (in.OpenSection("LongRange")){
    HasLongRange = true;
    ReadLongRangeHDF5(in);
    in.CloseSection();
  }
  else
    HasLongRange=false;

  // Print and debug
  //double U, dU, v;
  //double x(0.0), end(3.0);
  //while (x<=end) {
  //  calcUsqz(0,x,0,0,U,dU,v);
  //  cout << x << " " << U << " " << dU << " " << v << endl;
  //  x += 0.1;
  //}
  //dVec r(3), rp(3);
  //r = 0.;
  //rp = 0.;
  //rp(0) = 0.5;
  //double rmag = sqrt(dot(r,r));
  //double rpmag = sqrt(dot(rp,rp));
  //double s2 = dot(r-rp, r-rp);
  //double q = 0.5*(rmag + rpmag);
  //double z = (rmag - rpmag);
  //calcUsqz(sqrt(s2),q,z,0,U,dU,v);
  //cout << sqrt(s2) << " " << q << " " << z << " " << U << " " << dU << " " << v << endl;
}


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

