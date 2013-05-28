#ifndef NAN_H
#define NAN_H

inline 
bool myIsNAN (double x)
{
  union 
  {
    double d;
    unsigned long long int l;
  } val;
  val.d = x;
  return ((val.l == (unsigned long long int)0xfff8000000000000ULL) ||
	  (val.l == (unsigned long long int)0x7ff8000000000000ULL));
}

#endif
