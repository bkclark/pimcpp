#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <iostream>
#include <cassert>
#include <cstdlib>

typedef enum {PERIODIC, FIXED_FIRST, FIXED_SECOND, FLAT, NATURAL} BCType;

template<typename T>
class BoundaryCondition
{
private:
  BCType BC;
  T Val;
public:
  inline BCType GetType() { return BC; }
  T GetVal () {
    assert (BC == FIXED_FIRST || BC == FIXED_SECOND);
    return Val;
  }
  BoundaryCondition (BCType bctype)
  {
    if (bctype == FLAT) {
      BC = FIXED_FIRST;
      Val = T();
    }
    else if (bctype == NATURAL) {
      BC = FIXED_SECOND;
      Val = T();
    }
    else if (bctype == PERIODIC) {
      BC = PERIODIC;
    }
    else if (bctype == FIXED_FIRST) {
      std::cerr << "You must specify the derivative value"
		<< " for FIXED_FIRST BC.\n";
      abort();
    }
    else if (bctype == FIXED_SECOND) {
      std::cerr << "You must specify the FIXED_SECOND BC.\n";
      abort();
    }
    else {
      std::cerr << "Unknown boundary conditions.\n";
      abort();
    }
  }
  BoundaryCondition (BCType bctype, T val)
  {
    if ((bctype != FIXED_FIRST) && (bctype != FIXED_SECOND)) {
      std::cerr << "Cannot fix derivatives with periodic, flat, " 
		<< "or natural boundary conditions.\n";
      abort();
    }
    Val = val;
    BC = bctype;
  }

};


#endif
