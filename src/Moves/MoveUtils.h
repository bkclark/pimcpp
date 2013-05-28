#ifndef MOVE_UTILITIES_H
#define MOVE_UTILITIES_H


#include "../Blitz.h"
// typedef TinyVector<scalar,NDIM> dVec;

// A number of functions used by several 
// Move and Action classes
// -jg

dVec crossprod(dVec v1, dVec v2);
dVec cross(dVec v1, dVec v2);

// return the dot product of two vectors
double dotprod(dVec vec1, dVec vec2);

// return the dot product of two vectors,
// normalizing by a given value
// i.e. gives the dot product of two unit vectors
double dotprod(dVec vec1, dVec vec2, double mag);

// returns the magnitude of a vector
double Mag(dVec v);

// normalizes a vector
// i.e. generates a unit vector
dVec Normalize(dVec v);

// Scales a vector by a given scalar
// now deprecated
dVec Scale(dVec v, double scale);

// **Re-Scales** a vector by a given scalar
// used to be called Scale()
dVec Renormalize(dVec v, double scale);

// Takes a vector u and a unit vector R
// and returns the vector aligned that is
// the component of u parallel to R, and
// perp, the component of u perpendicular to R
void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);

// not really supported; just a quick hack
double GetAngle(dVec v1, dVec v2);
dVec GetBisector(dVec v1, dVec v2);

#endif
