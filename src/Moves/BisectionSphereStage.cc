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

#include <Common/MPI/Communication.h>
#include "BisectionSphereStage.h"

void BisectionStageSphereClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("SphereRadius",SphereRadius));


}

void BisectionStageSphereClass::CartesianToSpherical(dVec &r,double &theta,double &phi)
{
  //double a=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  theta = abs(atan2(sqrt(r[0]*r[0]+r[1]*r[1]),r[2]));
  phi = atan2(r[1],r[0]);
  //cerr <<"theta and phi"<<theta<<" "<<phi<<endl;
  //cerr<<"...x= "<<r[0]<<" y="<<r[1]<<" z="<<r[2]<<endl;
  return;
}

TinyVector<double,3> BisectionStageSphereClass::SphericalToCartesian(double &a,double &theta, double &phi)
{
  TinyVector<double,3> r;
  double sinTheta;
  double cosTheta;
  double sinPhi;
  double cosPhi;
  cosTheta=cos(theta);
  sinTheta=sin(theta);
  //  sincos(theta,&sinTheta,&cosTheta);
  //  sincos(phi,&sinPhi,&cosPhi);
  cosPhi=cos(phi);
  sinPhi=sin(phi);
  r[0]=a * sinTheta * cosPhi;
  r[1]=a * sinTheta * sinPhi;
  r[2]=a * cosTheta;
  return r;
}


void BisectionStageSphereClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  OutSection.FlushFile();
}

void BisectionStageSphereClass::Accept()
{
  //do nothing for now
  
}

void BisectionStageSphereClass::Reject()
{
  //do nothing for now

}

void BisectionStageSphereClass::SurfaceUnitVectors(dVec &Er,dVec &Etheta, dVec &Ephi, double theta, double phi)
{
  double sinPhi;
  double cosPhi;
  //  sincos(phi,&sinPhi,&cosPhi);
  cosPhi=cos(phi);
  sinPhi=sin(phi);
  double sinTheta;
  double cosTheta;
  //  sincos(theta,&sinTheta,&cosTheta);
  cosTheta=cos(theta);
  sinTheta=sin(theta);
  Er[0]=sinTheta*cosPhi;
  Er[1]=sinTheta*sinPhi;
  Er[2]=cosTheta;
  Etheta[0]=cosTheta*cosPhi;
  Etheta[1]=cosTheta*sinPhi;
  Etheta[2]=-sinTheta;

  Ephi[0]=-sinPhi;
  Ephi[1]=cosPhi;
  Ephi[2]=0;
}

///Projects the vector r onto the sphere of radius a
void BisectionStageSphereClass::ProjectOntoSphere(dVec &r,double a)
{
  double magR=sqrt(dot(r,r));
  dVec unitr;
  for (int dim=0;dim<NDIM;dim++){
    unitr(dim)=r(dim)/magR;
  }
  r=unitr*a;
  assert(dot(r,r)-a*a<1e-5);

}
double BisectionStageSphereClass::angleInBetween(dVec &r1,dVec &r2)
{
  double magr1=sqrt(dot(r1,r1));
  double magr2=sqrt(dot(r2,r2));
  double angle= acos(dot(r1,r2)/(magr1*magr2));
  //  if (angle<0)
  //    angle=angle+2*M_PI;
  return abs(angle);

}


double BisectionStageSphereClass::Sample(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  
  PathClass &Path = PathData.Path;
  int skip = 1<<(BisectionLevel+1);
  double levelTau = 0.5*PathData.Path.tau*skip;
  double randomTheta;
  double oldAngle;
  int numImages = PathData.Actions.NumImages;

  double logSampleProb=0.0;
  double logOldSampleProb=0.0;

  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    double a=SphereRadius;
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(2.0*lambda*levelTau/(a*a));
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
   
    for (int slice=slice1;slice<slice2;slice+=skip){
      SetMode(OLDMODE);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      ProjectOntoSphere(rbarOld,a);
      dVec rppOld=Path(slice+(skip>>1),ptcl);

      oldAngle=angleInBetween(rppOld,rbarOld);
      double randomDist;

      double thetaBar,phiBar;
      dVec rpp,rbar;
      SetMode(NEWMODE);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
      rbar=r+ 0.5*rdiff;

      ProjectOntoSphere(rbar,a);
      CartesianToSpherical(rbar,thetaBar,phiBar);
      dVec Er, Etheta, Ephi;
      SurfaceUnitVectors(Er,Etheta,Ephi,thetaBar,phiBar);

      randomTheta=abs(Path.Random.LocalGaussian(sigma));
      randomDist=a*randomTheta;

      double rn1=-1+2*Path.Random.Local();
      double rn2=-1+2*Path.Random.Local();
      double cosRandTheta,sinRandTheta;
      cosRandTheta=cos(randomTheta);
      sinRandTheta=sin(randomTheta);
      //      sincos(randomTheta,&sinRandTheta,&cosRandTheta);
      dVec samVec=rn1*Etheta+rn2*Ephi;
      rpp=a*cosRandTheta*Er+a*sinRandTheta*(samVec)/sqrt(dot(samVec,samVec));

      if (abs(sqrt(dot(rpp,rpp)) -SphereRadius)>pow(10.0,-5)) 
	cerr <<"Inside sphere!! "<<sqrt(dot(rpp,rpp))<<endl;
      ProjectOntoSphere(rpp,a);  
      

     ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
 
     
//       logSampleProb += prefactorOfSampleProb + 
// 	-0.5*(randomTheta*randomTheta)/sigma2+log(randomTheta*a);
//       logOldSampleProb += prefactorOfSampleProb + 
// 	-0.5*(oldAngle*oldAngle)/sigma2+log(oldAngle*a);
//      cerr<<randomTheta<<" "<<oldAngle<<" "<<abs(oldAngle)<<" "<<log(abs(randomTheta))<<" "<<log(abs(oldAngle))<<endl;

      logSampleProb += prefactorOfSampleProb + 
	-0.5*(randomTheta*randomTheta)/sigma2;
	//-log(abs(randomTheta));
      if (oldAngle==0)
	return 100;
	  
      else
	logOldSampleProb += prefactorOfSampleProb + 
	  -0.5*(oldAngle*oldAngle)/sigma2;
	  //-log(abs(oldAngle));
      
    }
  }

  return abs(oldAngle)/abs(randomTheta)*exp(-logSampleProb+logOldSampleProb);

}
// //////
// /////
// ////not saads way

// void BisectionStageSphereClass::CartesianToSpherical(dVec &r,double &theta,double &phi)
// {
//   double a=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
//   theta=atan2(sqrt(r[1]*r[1]+r[0]*r[0]),r[2]);
//   if (r[0]==0 && r[1]==0)
//     phi=0;
//   else
//     phi=atan2(r[1],r[0]);
//   //  if (r[0]<0)
//   //    phi=M_PI+phi;
  
//   //   if (theta<0)
// //     theta=theta+M_PI;
// //   if (phi<0)
// //     phi=phi+2*M_PI;
//   return;
// }


// // void BisectionStageSphereClass::CartesianToSpherical2(dVec &r,double &theta,double &phi)
// // {
// //   double a=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
// //   theta=atan2(r[2],sqrt(r[1]*r[1]+r[0]*r[0]));
// //   if (r[0]==0 && r[1]==0)
// //     phi=0;
// //   else
// //     phi=atan2(r[1],r[0]);
// //   //  if (r[0]<0)
// //   //    phi=M_PI+phi;

// // //   if (theta<0)
// // //     theta=theta+M_PI;
// // //   if (phi<0)
// // //     phi=phi+2*M_PI;
// //   return;
// // }

// TinyVector<double,3> BisectionStageSphereClass::SphericalToCartesian(double &a,double &theta, double &phi)
// {
//   TinyVector<double,3> r;
//   double sinTheta;
//   double cosTheta;
//   double sinPhi;
//   double cosPhi;
//   sincos(theta,&sinTheta,&cosTheta);
//   sincos(phi,&sinPhi,&cosPhi);
//   r[0]=a * sinTheta * cosPhi;
//   r[1]=a * sinTheta * sinPhi;
//   r[2]=a * cosTheta;
//   if (r[0]==-0)
//     r[0]=0;
//   if (r[1]==-0)
//     r[1]=0;
//   if (r[2]==-0)
//     r[2]=0;
	
//   return r;
// }


// ///Projects the vector r onto the sphere of radius a
// void BisectionStageSphereClass::ProjectOntoSphere(dVec &r,double a)
// {
//   double magR=sqrt(dot(r,r));
//   dVec unitr;
//   for (int dim=0;dim<NDIM;dim++){
//     unitr(dim)=r(dim)/magR;
//   }
//   r=unitr*a;
//   assert(dot(r,r)-a*a<1e-5);

// }

// void BisectionStageSphereClass::WriteRatio()
// { 
//   AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
//   OutSection.FlushFile();
// }

// void BisectionStageSphereClass::Accept()
// {
//   //do nothing for now
  
// }

// void BisectionStageSphereClass::Reject()
// {
//   //do nothing for now

// }
// double BisectionStageSphereClass::angleInBetween(dVec &r1,dVec &r2)
// {
//   double magr1=sqrt(dot(r1,r1));
//   double magr2=sqrt(dot(r2,r2));
//   double angle= acos(dot(r1,r2)/(magr1*magr2));
//   //  if (angle<0)
//   //    angle=angle+2*M_PI;
//   return abs(angle);

// }
// void BisectionStageSphereClass::RotateAroundZ(dVec &vec, double phi)
// {
//   phi=-phi;
//   if (phi<0)
//     phi=phi+2*M_PI;
//   double sinPhi;
//   double cosPhi;
//   sincos(phi,&sinPhi,&cosPhi);
//   double x=vec[0];
//   double y=vec[1];
//   double z=vec[2];
//   vec[0]=cosPhi * x + sinPhi * y;
//   vec[1]=-sinPhi*x+cosPhi*y;
//   vec[2]=z;
// }

// void BisectionStageSphereClass::RotateAroundX(dVec &vec, double theta)
// {
//   //  if (theta<0)
//   //    theta=theta+2*M_PI;
//   double sinTheta;
//   double cosTheta;
//   sincos(theta,&sinTheta,&cosTheta);
//   double x=vec[0];
//   double y=vec[1];
//   double z=vec[2];
//   vec[0]=x;
//   vec[1]=cosTheta * y+ sinTheta * z;
//   vec[2]=-(sinTheta) * y + cosTheta * z;
// }

// void BisectionStageSphereClass::RotateAroundY(dVec &vec, double theta)
// {
//   double sinTheta;
//   double cosTheta;
//   sincos(theta,&sinTheta,&cosTheta);
//   double x=vec[0];
//   double y=vec[1];
//   double z=vec[2];

//   vec[0]=cosTheta * x+ sinTheta * z;
//   vec[1]=y;
//   vec[2]=-sinTheta * x + cosTheta * z;
// }

// double BisectionStageSphereClass::Sample(int &slice1,int &slice2,
// 				   Array<int,1> &activeParticles)
// {
//   //  cerr<<"I'm being asked to sample"<<endl;
//   PathClass &Path = PathData.Path;
//   int skip = 1<<(BisectionLevel+1);
//   double levelTau = 0.5*PathData.Path.tau*skip;

//   int numImages = PathData.Actions.NumImages;

//   double logSampleProb=0.0;
//   double logOldSampleProb=0.0;

//   for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
//     int ptcl=activeParticles(ptclIndex);
//     double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
//     double sigma2=(1.0*lambda*levelTau);
//     double sigma=sqrt(sigma2);
//     double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
//     double a=SphereRadius;
//     for (int slice=slice1;slice<slice2;slice+=skip){
//       SetMode(OLDMODE);
//       dVec rOld=Path(slice,ptcl);
//       dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
//       dVec rbarOld=rOld+ 0.5*rdiffOld;

 
//       ProjectOntoSphere(rbarOld,a);
//       dVec rppOld=Path(slice+(skip>>1),ptcl);
//       double oldDist;
//       oldDist=angleInBetween(rppOld,rbarOld)*a;
//       //      if (BisectionLevel==1)
//       //	cerr<<"My old dist is "<<oldDist<<endl;
//       double randomDist;
//       double randomTheta;
//       double randomPhi;
//       double thetaBar,phiBar;
//       dVec rpp,rbar;
//       SetMode(NEWMODE);
//       dVec r=Path(slice,ptcl);
//       dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
//       rbar=r+ 0.5*rdiff;
      
   
//       ProjectOntoSphere(rbar,a);
//       CartesianToSpherical(rbar,thetaBar,phiBar);
//       randomDist=Path.Random.LocalGaussian(sigma);
//       Vec2 guassianVec;
//       Path.Random.LocalGaussianVec(sigma,guassianVec);
//       randomDist=sqrt(dot(guassianVec,guassianVec));
//       randomTheta=abs(randomDist/a);
//       randomPhi=Path.Random.Local()*2*M_PI;
//       //      if (BisectionLevel==1)
//       //	cerr<<"My random dist is "<<randomDist<<endl;
      
//       rpp=SphericalToCartesian(a,randomTheta,randomPhi);
//       RotateAroundY(rpp, thetaBar);
//       RotateAroundZ(rpp, phiBar);


      
      
// //       //Test code

// //       double modRandomTheta;
// //       modRandomTheta=randomTheta;
// //       //      if (abs(angleInBetween(rpp,rbar)-modRandomTheta)>1e-5){
// //       cerr<<"Angle is between: "<<angleInBetween(rpp,rbar)<<" "<<modRandomTheta<<endl;	
// //       dVec rbarDup;
// //       rbarDup[0]=0;
// //       rbarDup[1]=0;
// //       rbarDup[2]=a;
// //       RotateAroundY(rbarDup,thetaBar);
// //       cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
// //       RotateAroundZ(rbarDup,phiBar);
// //       //      cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
// //       cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
// //       cerr<<"Random dist:  "<<randomDist<<endl;
// //       assert(abs(angleInBetween(rpp,rbar)-modRandomTheta)<1e-5);
// //       //      }

// //       ///End Test Code
      
//       ////////////
//      ///Here we've stored the new position in the path
//       Path.SetPos(slice+(skip>>1),ptcl,rpp);
 
 
     
//       logSampleProb += prefactorOfSampleProb + 
// 	-0.5*(randomDist*randomDist)/sigma2;
//       logOldSampleProb += prefactorOfSampleProb + 
// 	-0.5*(oldDist*oldDist)/sigma2;

      
//     }
//   }

//   return exp(-logSampleProb+logOldSampleProb);

// }

