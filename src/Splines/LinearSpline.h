#ifndef LINEAR_SPLINE_H
#define LINEAR_SPLINE_H
#include <blitz/array.h>

class LinearSpline
{
private:
  double GridStart, GridEnd;
  blitz::Array<double,1> Data;
  double dx, dxInv;
public:
  inline double Start() { return GridStart; }
  inline double End()   { return GridEnd; }
  inline void Init (double start, double end, Array<double,1> &data);
  inline double operator()(double x);
};


inline double
LinearSpline::operator()(double x)
{
  double id = dxInv*(x-GridStart);
  double ipart;
  double frac = modf (id, &ipart);
  int i = (int) ipart;
  return ((1.0-frac)*Data(i) + frac*Data(i+1));
}

inline void
LinearSpline::Init (double start, double end, Array<double,1> &data)
{
  Data.resize(data.size());
  Data = data;
  GridStart = start;
  GridEnd   = end;
  dx = (GridEnd-GridStart)/(double)(data.size()-1);
  dxInv = 1.0/dx;
}
  


#endif
