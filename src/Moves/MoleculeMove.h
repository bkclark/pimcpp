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

#ifndef MOLECULE_MOVE_H
#define MOLECULE_MOVE_H

#include "MoleculeMoveBase.h"
#include "../Moves/MoveUtils.h"

class MoleculeRotate : public MolMoveClass
{
  bool doAllSlices;
  int sizeOfMulti;
 public:
  double Theta;
  void Set(double setTheta);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  MoleculeRotate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			cerr << "MoleculeRotate constructor" << endl;
    }

  //  double AcceptanceRatio(int numAccepted,int numMoves);

  //void WriteRatio()
  //  {
  //  };
};

class PairActionTestClass : public MolMoveClass
{
 public:
  double r0, rc;
  string filename;
  string spec1, spec2;
  int ptcl1, ptcl2;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  PairActionTestClass(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
    }

};

class BondStretch : public MolMoveClass
{
  bool doAllSlices;
 public:
  double s;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  BondStretch(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			cerr << "BondStretch constructor" << endl;
    }

  //  double AcceptanceRatio(int numAccepted,int numMoves);

  //void WriteRatio()
  //  {
  //  };
};

class MoleculeTranslate : public MolMoveClass
{
  bool doAllSlices;
	int counter;
 	public:
  double Step;
  void Set(double setStep);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  //inline void WriteRatio()
  //  {
  //  };

  MoleculeTranslate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }

};

// individual intramolecular moves
class ParticleTranslate : public MolMoveClass
{
	int counter;
 	public:
  double Sigma;
  void Set(double setSigma);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  //inline void WriteRatio()
  //  {
  //  };

  ParticleTranslate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }
};

class DimerMove : public MolMoveClass
{
 	public:
  double Step;
  void Set(double setStep);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  DimerMove(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
    }

};

class DummyEvaluate : public MolMoveClass
{
	int counter;
 	public:
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //inline void WriteRatio()
  //  {
  //  };

  DummyEvaluate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }
};

class MoleculeMulti : public MolMoveClass
{
  MoleculeRotate Rotate;
  MoleculeTranslate Trans;
  BondStretch Stretch;

 public:
  double Theta;
  double Step;
  double s;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  MoleculeMulti(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection),
    Rotate(myPathData, outSection), Trans(myPathData, outSection), Stretch(myPathData, outSection)
    {
			cerr << "MoleculeMulti constructor" << endl;
    }

};

#endif
