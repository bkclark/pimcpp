#ifndef SPACE_WARP_H
#define SPACE_WARP_H

#include <Common/Blitz.h>

class SpaceWarpClass
{
protected:
  int N;
  /// The current ion positions
  Array<Vec3,1> Rions;
  /// The proposed change in ion positions
  Array<Vec3,1> DeltaRions;
  Array<Vec3,1> rhat, Wgr;
  Array<double,1> Weights, g;
  Vec3 Box, BoxInv;
  bool DoForward;
  void PutInBox(Vec3 &r);
  /// The unnormalized weight function
  inline double phi(double rinv)      { return rinv*rinv*rinv*rinv;}
  /// Its logarithmic derivative
  inline double d_ln_phi(double rinv) { return -4.0*rinv; };
public:
  /// Set the ions and their change
  void Set (const Array<Vec3,1> &rions, const Array<Vec3,1> &delta, Vec3 box,
	    bool doForward=true);
  Vec3 ForwardWarp(Vec3 r, Mat3 &jmat);
  Vec3 ForwardWarp(Vec3 r);
  /// For a given ending position, this function returns where I must
  /// have warped from.  Solves iteratively.
  Vec3 ReverseWarp(Vec3 rp, Mat3 &jmat);
  Vec3 ReverseWarp(Vec3 rp);
  /// Warp based in the current set direction.
  inline Vec3 Warp (Vec3 r, Mat3 &jmat);
  /// This function sets r1p to make {r0,r1,r2} and {r0p,r1p,r2p}
  /// similar triangles.  It returns the length ratio |r2p-r0p|/|r2-r0|.
  void SimilarTriangles (const Vec3 &r0,  const Vec3 &r1,  const Vec3 &r2,
			 const Vec3 &r0p,       Vec3 &r1p, const Vec3 &r2p,
			 double &alpha, double &beta, double &h, 
			 double &ratio);
  /// Adjusts r1 to scale the height by factor s.  Returns the scaled height.
  double ScaleTriangleHeight (const Vec3 &r0, Vec3 &r1, const Vec3 &r2,
			    double s);
  // Solves the equation As^2 + B ln(s) + C = 0;  returns s.
  double SolveScaleEquation (double A, double B, double C);
};

inline Vec3
SpaceWarpClass::Warp (Vec3 r, Mat3 &jmat)
{
  if (DoForward)
    return ForwardWarp (r, jmat);
  else
    return ReverseWarp (r, jmat);
}


#endif
