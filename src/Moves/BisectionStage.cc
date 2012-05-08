/////////////////////////////////////////////////////////////
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


#include "BisectionStage.h"
#include "../MatrixOps/MatrixOps.h"


// double BisectionStageClass::LogSampleProb(int slice1, int slice2, 
// 					  Array<int,1> particles)
// {
//   PathClass &Path = PathData.Path;
//   dVec rpp;
//   int skip = 1<<(BisectionLevel+1);
//   double logSampleProb=0.0;

//   double levelTau = 0.5*PathData.Action.tau*skip;
//   for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
//     int ptcl = particles(ptclIndex);
//     double lambda=Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
//     for (int slice=slice1;slice<slice2;slice+=skip){
//       dVec r = Path(slice,ptcl);
//       dVec rdiff = Path.Velocity(slice, slice+skip, ptcl);

//       dVec rpp=Path(slice+(skip>>1),ptcl);
      
//       // We've ignored boundary conditions here (well we think this is
//       // fixed but we're not sure)  
//       dVec rbar=r + 0.5*rdiff;
//       dVec Delta= rpp - rbar;
//       Path.PutInBox(Delta);
      
//       double GaussProd=1.0;
//       for (int dim=0; dim<NDIM; dim++) {
// 	double GaussSum = 0.0;
// 	int NumImage = 1;
// 	for (int image=-NumImage; image <= NumImage; image++) {
// 	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
// 	  GaussSum += exp(-0.5*dist*dist/sigma2);
// 	}
// 	GaussProd *= GaussSum;
//       }
//       logSampleProb += prefactorOfSampleProb + log(GaussProd);

//     }
//   }
//   return logSampleProb;
// }

// double BisectionStageClass::SamplePaths(int startSlice, int endSlice, 
// 					Array<int,1> particles) 
// {
//   dVec rpp;
//   int skip = 1<<(BisectionLevel+1);
//   double logNewSampleProb=0.0;
//   PathClass &Path = PathData.Path;
//   ActionClass &Action = PathData.Action;
//   double levelTau = 0.5*PathData.Action.tau*skip;
//   for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
//     int  ptcl=particles(ptclIndex);
//     double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
//     for (int slice=startSlice;slice<endSlice;slice+=skip){
//       dVec r = Path(slice,ptcl);
//       //\\     dVec rp= Path(slice+skip,ptcl);
//       dVec rpp; //\\=Path(slice+(skip>>1),ptcl);
//       //      Path.PutInBox(r);
//       //      Path.PutInBox(rp);
//       //      Path.PutInBox(rpp);
//       dVec rdiff=Path.Velocity(slice, slice+skip, ptcl);
//       dVec rbar = r + 0.5*rdiff;
//       dVec newDelta;
//       Path.Random.LocalGaussianVec(sigma,newDelta);
//       PathData.Path.PutInBox(newDelta);
//       double GaussProd=1.0;
//       for (int dim=0; dim<NDIM; dim++) {
// 	  double GaussSum = 0.0;
// 	  int NumImage = 1;
// 	  for (int image=-NumImage; image <= NumImage; image++) {
// 	    double dist = newDelta[dim]+(double)image*Path.GetBox()[dim];
// 	    GaussSum += exp(-0.5*dist*dist/sigma2);
// 	  }
// 	  GaussProd *= GaussSum;
//       }
//       logNewSampleProb += prefactorOfSampleProb + log(GaussProd);
//       rpp=rbar+newDelta;

//       ///Here we've stored the new position in the path
//       Path.SetPos(slice+(skip>>1),ptcl,rpp);
//     }
//   } 
//   //  cerr<<"My logNewSampleProb is "<<logNewSampleProb<<endl;
//   return logNewSampleProb;
// }


//void BisectionStageClass::Read(IOSectionClass& IO)
//{




//}


void BisectionStageClass::CalcShift(Array<int,1> &activeParticles,int slice,
				    double sigma,
				    double &det)
{
  ///Correlated sampling
  cerr<<"Entering array"<<endl;
  int numActivePtcl=activeParticles.size();
  cerr<<"Active particle size is "<<numActivePtcl<<endl;
  Correlated.resize(NDIM*numActivePtcl,NDIM*numActivePtcl);
  Correlated=0.0;
  
  //  Array<dVec,1> dispShift(numActivePtcl);
  dispShift.resize(numActivePtcl);
  dispShift=0.0;
  //  double sigma;
  for (int ptclIndex=0;ptclIndex<numActivePtcl;ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double t1sum=0.0;
    for (int cptcl=0;cptcl<PathData.Path.NumParticles();cptcl++){
      double dist;
      dVec disp;
      PathData.Path.DistDisp(slice,ptcl,cptcl,dist,disp);
      double t1=PathData.Actions.ShortRange.dUdR(slice,ptcl,cptcl,BisectionLevel);
      double t2=PathData.Actions.ShortRange.d2UdR2(slice,ptcl,cptcl,BisectionLevel);
      for (int dim=0;dim<NDIM;dim++){
	dispShift(ptclIndex)[dim]-=t1*disp[dim];
	Correlated(ptclIndex*NDIM+dim,ptclIndex*NDIM+dim)-=t2*disp[dim]*disp[dim];
	for (int dim2=0;dim2<dim;dim2++){
	  Correlated(ptclIndex*NDIM+dim,ptclIndex*NDIM+dim2)-=t2*disp[dim]*disp[dim2];
	  Correlated(ptclIndex*NDIM+dim2,ptclIndex*NDIM+dim)-=t2*disp[dim]*disp[dim2];
	}
      }
      t1sum+=t1;
    }
    t1sum*=sigma;
    for (int dim=0;dim<NDIM;dim++)
      Correlated(ptclIndex*NDIM+dim,ptclIndex*NDIM+dim)+=t1sum;
    for (int ptclIndex2=0;ptclIndex2<numActivePtcl;ptclIndex2++){
      int ptcl2=activeParticles(ptclIndex2);
      double dist;
      dVec disp;
      PathData.Path.DistDisp(slice,ptcl,ptcl2,dist,disp);
      double t2=PathData.Actions.ShortRange.d2UdR2(slice,ptcl,ptcl2,BisectionLevel);
      for (int dim=0;dim<NDIM;dim++){
	Correlated(ptclIndex*NDIM+dim,ptclIndex2*NDIM+dim)-=t2*disp[dim]*disp[dim];
	for (int dim2=0;dim2<NDIM;dim2++){
	  Correlated(ptclIndex*NDIM+dim,ptclIndex2*NDIM+dim2)-=t2*disp[dim]*disp[dim2];
	  Correlated(ptclIndex*NDIM+dim2,ptclIndex2*NDIM+dim)-=t2*disp[dim]*disp[dim2];
	  Correlated(ptclIndex*NDIM+dim2,ptclIndex2*NDIM+dim)-=t2*disp[dim]*disp[dim2];
	  Correlated(ptclIndex*NDIM+dim,ptclIndex2*NDIM+dim2)-=t2*disp[dim]*disp[dim2];
	}
	  
      }
    }
  }
  cerr<<"Printing dispShift"<<endl;
  for (int i=0;i<dispShift.size();i++)
    cerr<<dispShift(i)<<endl;
  cerr<<"corrlated array"<<endl;
  for (int i=0;i<Correlated.extent(0);i++){
    for (int j=0;j<Correlated.extent(1);j++)
      cerr<<Correlated(i,j)<<" ";
    cerr<<endl;
  }
  cerr<<"Done with array"<<endl;
  S.resize(Correlated.extent(0),Correlated.extent(1));
  for (int i=0;i<Correlated.extent(0);i++)
    for (int j=0;j<Correlated.extent(1);j++)
      S(i,j)=Correlated(i,j);
  det=Determinant (S);
  //CholeskyBig(S);
  //cerr<<"corrlated Cholesky array"<<endl;
  cerr<<"ERROR no call to Cholesky array"<<endl;
  assert(0);
  for (int i=0;i<S.extent(0);i++){
    for (int j=0;j<S.extent(1);j++)
      cerr<<S(i,j)<<" ";
    cerr<<endl;
  }
  

}


void BisectionStageClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

void BisectionStageClass::Accept()
{
  LocalStageClass::Accept();
  //do nothing for now
  
}

void BisectionStageClass::Reject()
{
  LocalStageClass::Reject();
  //do nothing for now

}

///Calculates a new rbar that is displaced
///by the old rbar via the correlated sampling
void guassianDisplace(dVec rbar)
{
  


}

double BisectionStageClass::Sample(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);


  PathClass &Path = PathData.Path;



  if (UseCorrelatedSampling){
    Correlated.resize(NDIM*activeParticles.size(),NDIM*activeParticles.size());
  }

  int skip = 1<<(BisectionLevel+1);
  double levelTau = 0.5*PathData.Path.tau*skip;

  int numImages = PathData.Actions.NumImages;
  int oldSlice2=slice2;
  slice2=slice1+(1<<TotalLevels);
  double logSampleProb=0.0;
  double logOldSampleProb=0.0;
  dVec firstDelta;

  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=slice1;slice<slice2;slice+=skip){
      SetMode(OLDMODE);
      if (UseCorrelatedSampling){
	      ///testing the calcShift 
	      //done testing
	      ///Correlated sampling
	      dVec rshiftOld;
	      rshiftOld=0.0;
	      for (int cptcl=0;cptcl<PathData.Path.NumParticles();cptcl++){
	        double dist;
	        dVec disp;
	        PathData.Path.DistDisp(slice,ptcl,cptcl,dist,disp);
	        double t1=PathData.Actions.ShortRange.dUdR(slice,ptcl,cptcl,BisectionLevel);
      
	        for (int dim=0;dim<NDIM;dim++)
	          rshiftOld[dim]-=t1*disp[dim];
	      }    
	      //	cerr<<"R shift old"<<endl;
	      //	cerr<<rshiftOld<<endl;
       
      }
      ///      cerr<<"Ptcl is "<<ptcl<<" "<<rshiftOld[0]<<" "<<rshiftOld[1]<<endl;
      ///


      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      ///correlated sampling
      if (UseCorrelatedSampling){
	      dVec rshiftOld;
	      double drmax=0.1*PathData.Path.GetBox()[0];
	      for (int dim=0;dim<NDIM;dim++)
	        if (-drmax<=rshiftOld[dim] && rshiftOld[dim]<=drmax)
	          rbarOld[dim] +=rshiftOld[dim];
	        else if (rshiftOld[dim] < -drmax)
	          rbarOld[dim] +=-drmax;
	        else 
	          rbarOld[dim] +=drmax;
      }
      ///
      dVec rppOld=Path(slice+(skip>>1),ptcl);
      dVec DeltaOld=rppOld-rbarOld;
      Path.PutInBox(DeltaOld);
      
      SetMode(NEWMODE);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
  
      dVec rbar=r+ 0.5*rdiff;

      if (UseCorrelatedSampling){
	      ///Correlated sampling
	      dVec rshift;
	      rshift=0.0;
	      for (int cptcl=0;cptcl<PathData.Path.NumParticles();cptcl++){
	        double dist;
	        dVec disp;
	        PathData.Path.DistDisp(slice,ptcl,cptcl,dist,disp);
	        double t1=PathData.Actions.ShortRange.dUdR(slice,ptcl,cptcl,BisectionLevel);
	        for (int dim=0;dim<NDIM;dim++)
	          rshift[dim]-=t1*disp[dim];
	      }     
	
	      double drmax=0.1*PathData.Path.GetBox()[0];
	      for (int dim=0;dim<NDIM;dim++)
	        if (-drmax<=rshift[dim] && rshift[dim]<=drmax)
	          rbar[dim] +=rshift[dim];
	        else if (rshift[dim] < -drmax)
	          rbar[dim] +=-drmax;
	        else 
	          rbar[dim] +=drmax;
      }

      ///Here we've stored the new position in the path
      dVec  Delta;
      Path.Random.LocalGaussianVec(sigma,Delta);
      PathData.Path.PutInBox(Delta);
//       if (ptclIndex==0)
//       	firstDelta=Delta;
//       if (ptclIndex==1){
//       	double length=sqrt(dot(Delta,Delta));
//       	double firstLength=sqrt(dot(firstDelta,firstDelta));
//       	Delta=-1.0*firstDelta*length/firstLength;
//       }
      double GaussProd=1.0;
      double GaussProdOld=1.0;
      for (int dim=0; dim<NDIM; dim++) {
	double GaussSum = 0.0;
	double GaussSumOld =0.0;
	for (int image=-numImages; image <= numImages; image++) {
	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
	  double distOld=DeltaOld[dim]+(double)image*Path.GetBox()[dim];
	  GaussSum += exp(-0.5*dist*dist/sigma2);
	  GaussSumOld += exp(-0.5*distOld*distOld/sigma2);
	}
	GaussProd *= GaussSum;
	GaussProdOld *= GaussSumOld;
      }
      logOldSampleProb += prefactorOfSampleProb + log(GaussProdOld);
      logSampleProb += prefactorOfSampleProb  + log(GaussProd);
      
      dVec rpp=rbar+Delta;
      ////REALLY BAD HACK!
      ///      dVec rpp=rbar+100*Delta;
      ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
    }
  }
//   SetMode (NEWMODE);
//   SamplePaths (slice1, slice2, activeParticles);
//   SetMode(OLDMODE);
//   double oldSample = LogSampleProb (slice1, slice2, activeParticles);
//   SetMode(NEWMODE);
//   double newSample = LogSampleProb (slice1, slice2, activeParticles);
  //cerr << "oldSample = " << oldSample << endl;
  //cerr << "oldSampleProb = " << logOldSampleProb << endl;
  //cerr << "newSample = " << newSample << endl;
  //cerr << "newSampleProb = " << logSampleProb << endl;
  //  if (activeParticles.size()>1)
  //    cerr<<"Value is "<<BisectionLevel<<": "<<exp(-logSampleProb+logOldSampleProb)<<endl;
  slice2=oldSlice2;
  //  NewSample=exp(logSampleProb);
  //  OldSample=exp(logOldSampleProb);


	gettimeofday(&end,   &tz);

	//	//	TimeSpent += (double)(end.tv_sec-start.tv_sec) +
	//	//	  1.0e-6*(double)(end.tv_usec-start.tv_usec);
	//	cerr<<"What's time of day"<<TimeSpent<<endl;
  return exp(-logSampleProb+logOldSampleProb);
  //return (exp (-newSample + oldSample));
}

double BisectionStageClass::Sample_old(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;
  int skip = 1<<(BisectionLevel+1);
  double levelTau = 0.5*PathData.Path.tau*skip;
  int numImages = PathData.Actions.NumImages;
  int oldSlice2=slice2;
  slice2=slice1+(1<<TotalLevels);
  double logSampleProb=0.0;
  double logOldSampleProb=0.0;
  Array<double,1> sChiOld(NDIM*activeParticles.size());
  ///BUG: will only work if all of activeParticles are of the same species
  for (int slice=slice1;slice<slice2;slice+=skip){
    double lambda=PathData.Path.ParticleSpecies(activeParticles(0)).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    SetMode(OLDMODE);
    double oldDet;
    CalcShift(activeParticles,slice,sigma,oldDet);
    logOldSampleProb+=log(sqrt(oldDet));
    for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
      int ptcl=activeParticles(ptclIndex);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+0.5*rdiffOld;
      double drmax=0.1*PathData.Path.GetBox()[0];
      for (int dim=0;dim<NDIM;dim++)
	if (-drmax<=dispShift(ptclIndex)(dim) && dispShift(ptclIndex)(dim)<=drmax)
	  rbarOld(dim) +=dispShift(ptclIndex)(dim);
	else if (dispShift(ptclIndex)(dim) < -drmax)
	  rbarOld(dim) +=-drmax;
	else 
	  rbarOld(dim) +=drmax;
      

      for (int dim=0;dim<NDIM;dim++)
	sChiOld(NDIM*ptclIndex+dim)=rOld(dim)-rbarOld(dim);
    }
    Array<double,1> AsChi(NDIM*activeParticles.size());
    GJInverse(Correlated);
    MatVecProd (Correlated, sChiOld,AsChi);
    double oldProbTerm=InnerProduct(AsChi,sChiOld);
    logOldSampleProb -=0.5*oldProbTerm;
    
    SetMode(NEWMODE);
    double newDet;
    CalcShift(activeParticles,slice,sigma,newDet);
    logSampleProb += log(sqrt(newDet));
    Array<double,1> gaussianVec(NDIM*activeParticles.size());
    for (int i=0;i<NDIM*activeParticles.size();i++)
      gaussianVec(i)=Path.Random.LocalGaussian(1.0);
    Array<double,1> sChi(NDIM*activeParticles.size());
    MatVecProd(S,gaussianVec,sChi);
    GJInverse(Correlated);
    MatVecProd(Correlated,sChi,AsChi);
    double newProbTerm=InnerProduct(AsChi,sChi);
    logSampleProb-=0.5*newProbTerm;
    for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
      int ptcl=activeParticles(ptclIndex);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbar=r+ 0.5*rdiff;
      double drmax=0.1*PathData.Path.GetBox()[0];
      for (int dim=0;dim<NDIM;dim++)
	if (-drmax<=dispShift(ptclIndex)(dim) && dispShift(ptclIndex)(dim)<=drmax)
	  rbar(dim) +=dispShift(ptclIndex)(dim);
	else if (dispShift(ptclIndex)[dim] < -drmax)
	  rbar(dim) +=-drmax;
	else 
	  rbar(dim) +=drmax;
      dVec rpp;
      for (int dim=0;dim<NDIM;dim++)
	rpp(dim)=rbar(dim)+sChi(ptclIndex*NDIM+dim);
      ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
    }
  }
  slice2=oldSlice2;

  return exp(-logSampleProb+logOldSampleProb);
}

