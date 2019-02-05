
// tictoc.h

#ifdef CLOCK_HIGHRES
#define SAMPLED_CLOCK CLOCK_HIGHRES
#else
#define SAMPLED_CLOCK CLOCK_REALTIME
#endif

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
