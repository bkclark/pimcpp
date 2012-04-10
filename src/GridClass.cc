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

#include "PathClass.h"
#include "GridClass.h"


///affected cells are the cells that are within the cutoff and need to be used.  
///Then GridsArray contains a list of the particles in any given cell



double minAbs(double d1, double d2)
{
  if (abs(d1)<abs(d2))
    return d1;
  else 
    return d2;

}

///Grids affect is used 
bool CellMethodClass::GridsAffect(CellInfoClass &grid1,CellInfoClass &grid2)
{
  
  dVec diff1=grid1.left-grid2.right;
  dVec diff2=grid1.right-grid2.left;
  Path.PutInBox(diff1);
  Path.PutInBox(diff2);
  dVec minDisp;
  for (int dim=0;dim<NDIM;dim++)
    minDisp[dim]=minAbs(diff1[dim],diff2[dim]);
  return (sqrt(dot(minDisp,minDisp))<CutoffDistance);
}


#if NDIM==3
//used
void CellMethodClass::Init(dVec box,Array<int,1> numGrid)
{
  cerr<<"Starting my initialization"<<endl;
  NumGrid.resize(NDIM);
  for (int dim=0;dim<NDIM;dim++){
    NumGrid(dim)=numGrid(dim);
    cerr<<"Num grid "<<dim<<" is"<<numGrid(dim)<<endl;
  }
  GridsArray.resize(numGrid(0),numGrid(1),numGrid(2));
  //  GridsArray.resize(NumGrid);
  double xStart=-box[0]/2.0;
  double yStart=-box[1]/2.0;
  double zStart=-box[2]/2.0;
  double xSize=box[0]/numGrid(0);
  double ySize=box[1]/numGrid(1);
  double zSize=box[2]/numGrid(2);
  Xeffect=(int)(ceil(CutoffDistance/xSize));
  Yeffect=(int)(ceil(CutoffDistance/ySize));
  Zeffect=(int)(ceil(CutoffDistance/zSize));
  cerr<<"My effect is "<<Xeffect<<" "<<Yeffect<<" "<<Zeffect<<endl;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      for (int zCnt=0;zCnt<numGrid(2);zCnt++){
	//	cerr<<"We have "<<xCnt<<" "<<yCnt<<" "<<zCnt<<endl;
	GridsArray(xCnt,yCnt,zCnt).left[0]=xCnt*xSize+xStart;
	GridsArray(xCnt,yCnt,zCnt).right[0]=(xCnt+1)*xSize+xStart;
	GridsArray(xCnt,yCnt,zCnt).left[1]=yCnt*ySize+yStart;
	GridsArray(xCnt,yCnt,zCnt).right[1]=(yCnt+1)*ySize+yStart;
	GridsArray(xCnt,yCnt,zCnt).left[2]=zCnt*zSize+zStart;
	GridsArray(xCnt,yCnt,zCnt).right[2]=(zCnt+1)*zSize+zStart;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(0)=xCnt;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(1)=yCnt;
	GridsArray(xCnt,yCnt,zCnt).MyLoc(2)=zCnt;
	GridsArray(xCnt,yCnt,zCnt).Particles.resize(Path.NumTimeSlices());				      
	//	cerr<<"Done"<<endl;
      }
    }
  }
  int numAffected=0;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      for (int zCnt=0;zCnt<numGrid(2);zCnt++){
	if (GridsAffect(GridsArray(xCnt,yCnt,zCnt),GridsArray(0,0,0)))
	  numAffected++;
      }
    }
  }
  cerr<<"HELLO! I'M HERE!!!!"<<endl;
  AffectedCells.resize(numAffected);
  int onVal=0;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      for (int zCnt=0;zCnt<numGrid(2);zCnt++){
	if (GridsAffect(GridsArray(xCnt,yCnt,zCnt),GridsArray(0,0,0))){
	  AffectedCells(onVal)[0]=xCnt;
	  AffectedCells(onVal)[1]=yCnt;
	  AffectedCells(onVal)[2]=zCnt;
	  onVal++;
	  cerr<<"Affected: "<<xCnt<<" "<<yCnt<<" "<<zCnt<<endl;
	}
      }
    }
  }
	
  cerr<<"Ending my initialization"<<endl;
}
#endif


#if NDIM==2
//used
void CellMethodClass::Init(dVec box,Array<int,1> numGrid)
{
  cerr<<"Starting my initialization"<<endl;
  NumGrid.resize(NDIM);
  for (int dim=0;dim<NDIM;dim++){
    NumGrid(dim)=numGrid(dim);
    cerr<<"Num grid "<<dim<<" is"<<numGrid(dim)<<endl;
  }
  GridsArray.resize(numGrid(0),numGrid(1));
  double xStart=-box[0]/2.0;
  double yStart=-box[1]/2.0;
  double xSize=box[0]/numGrid(0);
  double ySize=box[1]/numGrid(1);
  Xeffect=(int)(ceil(CutoffDistance/xSize));
  Yeffect=(int)(ceil(CutoffDistance/ySize));
  cerr<<"My effect is "<<Xeffect<<" "<<Yeffect<<" "<<endl;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      GridsArray(xCnt,yCnt).left[0]=xCnt*xSize+xStart;
      GridsArray(xCnt,yCnt).right[0]=(xCnt+1)*xSize+xStart;
      GridsArray(xCnt,yCnt).left[1]=yCnt*ySize+yStart;
      GridsArray(xCnt,yCnt).right[1]=(yCnt+1)*ySize+yStart;
      GridsArray(xCnt,yCnt).MyLoc(0)=xCnt;
      GridsArray(xCnt,yCnt).MyLoc(1)=yCnt;
      GridsArray(xCnt,yCnt).Particles.resize(Path.NumTimeSlices());				      
    }
  }
  int numAffected=0;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
      if (GridsAffect(GridsArray(xCnt,yCnt),GridsArray(0,0)))
	numAffected++;
    }
  }
  cerr<<"HELLO! I'M HERE!!!!"<<endl;
  AffectedCells.resize(numAffected);
  int onVal=0;
  for (int xCnt=0;xCnt<numGrid(0);xCnt++){
    for (int yCnt=0;yCnt<numGrid(1);yCnt++){
	if (GridsAffect(GridsArray(xCnt,yCnt),GridsArray(0,0))){
	  AffectedCells(onVal)[0]=xCnt;
	  AffectedCells(onVal)[1]=yCnt;
	  onVal++;
	  cerr<<"Affected: "<<xCnt<<" "<<yCnt<<" "<<endl;
	}
    }
  }
  cerr<<"Ending my initialization"<<endl;
}
#endif


#if NDIM==3
///This needs to be rewritten to be a reasonable speed!!
void CellMethodClass::FindBox(dVec myPoint,int &x,int &y,int &z)
{
  Path.PutInBox(myPoint);
  dVec box=Path.GetBox();
  double xStart=-box[0]/2.0;
  double yStart=-box[1]/2.0;
  double zStart=-box[2]/2.0;
  double xSize=box[0]/NumGrid(0);
  double ySize=box[1]/NumGrid(1);
  double zSize=box[2]/NumGrid(2);
  x=(int)floor((myPoint[0]-xStart-0.001)/xSize);
  y=(int)floor((myPoint[1]-yStart-0.001)/ySize);
  z=(int)floor((myPoint[2]-zStart-0.001)/zSize);
 x=x+(x<0);
  y=y+(y<0);
  z=z+(z<0);

}
#endif  

#if NDIM==2
///This needs to be rewritten to be a reasonable speed!!
void CellMethodClass::FindBox(dVec myPoint,int &x,int &y)
{
  Path.PutInBox(myPoint);
  dVec box=Path.GetBox();
  dVec box_inv=Path.GetBoxInv();
  double xStart=-box[0]*0.5; 
  double yStart=-box[1]*0.5; // /2.0;
  double one_over_xSize=NumGrid(0)*box_inv(0);
  double one_over_ySize=NumGrid(1)*box_inv(1);
  x=(int)floor((myPoint[0]-xStart-0.001)*one_over_xSize);
  y=(int)floor((myPoint[1]-yStart-0.001)*one_over_ySize);
  x=x+(x<0);
  y=y+(y<0);

}
#endif  

#if NDIM==3
///used
void CellMethodClass::BinParticles(int slice)
{
  int x,y,z;
  for (int x=0;x<NumGrid(0);x++)
    for (int y=0;y<NumGrid(1);y++)
      for (int z=0;z<NumGrid(2);z++)
	GridsArray(x,y,z).Particles(slice).clear();
  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
    FindBox(Path(slice,ptcl),x,y,z);
    GridsArray(x,y,z).Particles(slice).push_back(ptcl);
  }	    
}
#endif


#if NDIM==2
///used
void CellMethodClass::BinParticles(int slice)
{
  int x,y,z;
  for (int x=0;x<NumGrid(0);x++)
    for (int y=0;y<NumGrid(1);y++)
      GridsArray(x,y).Particles(slice).clear();
  for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
    FindBox(Path(slice,ptcl),x,y);
    GridsArray(x,y).Particles(slice).push_back(ptcl);
  }	    
}
#endif

#if NDIM==3
void CellMethodClass::ReGrid(int slice,int ptcl)
{

  int currX,currY,currZ;
  int newX,newY,newZ;
  SetMode(OLDMODE);
  FindBox(Path(slice,ptcl),currX,currY,currZ);
  SetMode(NEWMODE);
  FindBox(Path(slice,ptcl),newX,newY,newZ);
  if (newX!=currX || newY!=currY || newZ!=currZ){
    GridsArray(currX,currY,currZ).Particles(slice).remove(ptcl);
    GridsArray(newX,newY,newZ).Particles(slice).push_back(ptcl);
  }
} 
#endif



#if NDIM==2
void CellMethodClass::ReGrid(int slice,int ptcl)
{

  int currX,currY;
  int newX,newY;
  SetMode(OLDMODE);
  FindBox(Path(slice,ptcl),currX,currY);
  SetMode(NEWMODE);
  FindBox(Path(slice,ptcl),newX,newY);
  if (newX!=currX || newY!=currY){
    GridsArray(currX,currY).Particles(slice).remove(ptcl);
    GridsArray(newX,newY).Particles(slice).push_back(ptcl);
  }
} 
#endif

#if NDIM==3
void CellMethodClass::PrintParticles(int slice)
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      for (int z=0;z<NumGrid(2);z++){
	cerr<<"I am grid: "<<x<<" "<<y<<" "<<z<<": ";
	for (list<int>::iterator i=GridsArray(x,y,z).Particles(slice).begin();i!=GridsArray(x,y,z).Particles(slice).end();i++){
	  cerr<<*i<<" ";
	}
	cerr<<endl;
      }
    }
  }
}	
#endif


#if NDIM==2
void CellMethodClass::PrintParticles(int slice)
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      cerr<<"I am grid: "<<x<<" "<<y<<" "<<": ";
      for (list<int>::iterator i=GridsArray(x,y).Particles(slice).begin();i!=GridsArray(x,y).Particles(slice).end();i++){
	cerr<<*i<<" ";
      }
      cerr<<endl;
    }
  }
}	
#endif

#if NDIM==3
void CellMethodClass::PrintNeighborGrids()
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      for (int z=0;z<NumGrid(2);z++){
	cerr<<"I am grid: "<<x<<" "<<y<<" "<<z;
	int dummyVar=5;
      }
    }
  }
}	
#endif								    

#if NDIM==2
void CellMethodClass::PrintNeighborGrids()
{
  for (int x=0;x<NumGrid(0);x++){
    for (int y=0;y<NumGrid(1);y++){
      cerr<<"I am grid: "<<x<<" "<<y;
      int dummyVar=5;
    }
  }
}	
#endif								    

			     
//deprecated
// bool CellMethodClass::InBox(CellInfoClass &theGrid,dVec thePoint) 
// {
//   Path.PutInBox(thePoint);
//   return (
//   	  (theGrid.left[0]<=thePoint[0] && thePoint[0]<theGrid.right[0]) &&
//   	  (theGrid.left[1]<=thePoint[1] && thePoint[1]<theGrid.right[1]) && 
//   	  (theGrid.left[2]<=thePoint[2] && thePoint[2]<theGrid.right[2])
//   	  );

// }




// ///Deprecated
// void CellMethodClass::BuildNeighborGrids()
// {
//   cerr<<"begin"<<endl;
//   cerr<<NumGrid(0)<<" "<<NumGrid(1)<<" "<<NumGrid(2)<<endl;
//   for (int x=0;x<NumGrid(0);x++){
//     for (int y=0;y<NumGrid(1);y++){
//       for (int z=0;z<NumGrid(2);z++){
// 	for (int x2=0;x2<NumGrid(0);x2++){
// 	  for (int y2=0;y2<NumGrid(1);y2++){
// 	    for (int z2=0;z2<NumGrid(2);z2++){
// 	      if (GridsAffect(GridsArray(x,y,z),GridsArray(x2,y2,z2))){
// 		int dummyVar=5;
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   cerr<<endl;
// }
