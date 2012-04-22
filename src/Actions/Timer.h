#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
class TimerClass
{
 public:
  struct timeval start, end;
  struct timezone tz;
  double TimeSpent;
  string Name;

 TimerClass(string tempName) : Name(tempName)
  {
    TimeSpent=0.;
    
  }


  void Start()
  {
    gettimeofday(&start, &tz);
  }

  void Stop()
  {
    gettimeofday(&end,   &tz);
    TimeSpent += (double)(end.tv_sec-start.tv_sec) +
      1.0e-6*(double)(end.tv_usec-start.tv_usec);
  }

  void Clear()
  {
    TimeSpent=0.0;
  }

  double Time()
  {
    return TimeSpent;

  }

};
#endif
