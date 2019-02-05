
// tictoc.cpp

#include "include.h"

#ifdef HAVE_CLOCK_GETTIME

void tic(tictoc* tv)
{
  if(clock_gettime(SAMPLED_CLOCK, tv) )
  {
    tv->tv_sec = tv->tv_nsec = -1;
  }
}

double toc(tictoc* tv)
{
  struct timespec tv2;

  memset(&tv2, 0, sizeof(struct timespec) );

  if(clock_gettime(SAMPLED_CLOCK, &tv2) )
  {
    tv2.tv_sec = tv2.tv_nsec = -1;
  }

  double sec = static_cast<double>(tv2.tv_sec - tv->tv_sec);

  double nsec = static_cast<double>(tv2.tv_nsec - tv->tv_nsec);

  return (sec + 1.0e-9 * nsec);
}

#elif HAVE_GETTIMEOFDAY

void tic(tictoc* tv)
{
  gettimeofday(tv, 0L);
}

double toc(tictoc* tv)
{
  tictoc tv2;

  memset(&tv2, 0, sizeof(tictoc) );

  gettimeofday(&tv2, 0L);

  double sec = static_cast<double>(tv2.tv_sec - tv->tv_sec);

  double usec = static_cast<double>(tv2.tv_usec - tv->tv_usec);

  return (sec + 1.0e-6 * usec);
}

#elif _WIN32

void tic(tictoc* tv)
{
  *tv = GetTickCount();
}

double toc(tictoc* tv)
{
  tictoc tv2;

  memset(&tv2, 0, sizeof(tictoc) );

  tv2 = GetTickCount();

  return 1.0e-3 * (tv2 - *tv);
}

#else

void tic(tictoc* tv)
{
  time(tv);
}

double toc(tictoc* tv)
{
  tictoc tv2;

  memset(&tv2, 0, sizeof(tictoc) );

  time(&tv2);

  return difftime(tv2, *tv);
}

#endif
