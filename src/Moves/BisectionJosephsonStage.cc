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


#include "BisectionJosephsonStage.h"

void BisectionJosephsonStageClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

void BisectionJosephsonStageClass::Accept()
{
  //do nothing for now
  
}

void BisectionJosephsonStageClass::Reject()
{
  //do nothing for now

}

double BisectionJosephsonStageClass::Sample(int &slice1,int &slice2,
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
  dVec firstDelta;
  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);

    for (int slice=slice1;slice<slice2;slice+=skip){

      

      SetMode(OLDMODE);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      dVec rppOld=Path(slice+(skip>>1),ptcl);
      dVec DeltaOld=rppOld-rbarOld;
      Path.PutInBox(DeltaOld);
      

      ///Sortof HACK
      SetMode(NEWMODE);
      if (PathData.Path.Random.Local()>0.9){
	if (BisectionLevel==2){
	  if (PathData.Path.Random.Local()>0.5){
	    Path(slice+(skip>>1),0)[0]=Path(slice+(skip>>1),0)[0]+2*M_PI;
	  }
	  else {
	    Path(slice+(skip>>1),0)[0]=Path(slice+(skip>>1),0)[0]-2*M_PI;
	  }
	}
	
      }


      ///End Sortof HACK



      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
  
      dVec rbar=r+ 0.5*rdiff;
      ///Here we've stored the new position in the path
      dVec  Delta;
      Path.Random.LocalGaussianVec(sigma,Delta);
      PathData.Path.PutInBox(Delta);
      Delta[1]=0; Delta[2]=0;
      double GaussProd=1.0;
      double GaussProdOld=1.0;
      //      HACK! for (int dim=0; dim<NDIM; dim++) {
//       for (int dim=0; dim<1; dim++) {
	double GaussSum = 0.0;
 	double GaussSumOld =0.0;
// 	for (int image=-numImages; image <= numImages; image++) {
// 	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
// 	  double distOld=DeltaOld[dim]+(double)image*Path.GetBox()[dim];
// 	  GaussSum += exp(-0.5*dist*dist/sigma2);
// 	  GaussSumOld += exp(-0.5*distOld*distOld/sigma2);
// 	}
// 	GaussProd *= GaussSum;
// 	GaussProdOld *= GaussSumOld;
//       }
//       logOldSampleProb += prefactorOfSampleProb + log(GaussProdOld);
//       logSampleProb += prefactorOfSampleProb  + log(GaussProd);

 
      double dist = Delta[0];
      double distOld=DeltaOld[0];
      GaussSum += (-0.5*dist*dist/sigma2);
      GaussSumOld += (-0.5*distOld*distOld/sigma2);
      logOldSampleProb =GaussSumOld;
      logSampleProb += GaussSum;
      //      cerr<<sigma2<<" "<<-0.5*dist*dist/sigma2;
      dVec rpp=rbar+Delta;
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
      //      cerr<<"My delta is "<<Delta<<endl;
    }
  }
  slice2=oldSlice2;
  //  cerr<<"I have sample probs of "<<-logSampleProb<<" "
  //      <<logOldSampleProb<<" "<<(-logSampleProb+logOldSampleProb)<<endl;
  return exp(-logSampleProb+logOldSampleProb);
  
}

