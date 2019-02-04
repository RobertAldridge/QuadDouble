
// tictoc.h

#ifdef HAVE_CLOCK_GETTIME
typedef struct timespec tictoc;
#elif HAVE_GETTIMEOFDAY
typedef struct timeval tictoc;
#elif _WIN32
typedef DWORD tictoc;
#else
typedef time_t tictoc;
#endif

void tic(tictoc* tv);

double toc(tictoc* tv);
