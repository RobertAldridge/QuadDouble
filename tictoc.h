
// tictoc.h

/*
 * Contains function used for timing.
 */

#ifdef HAVE_CLOCK_GETTIME
typedef struct timespec tictoc;
#elif HAVE_GETTIMEOFDAY
typedef struct timeval tictoc;
#elif _WIN32
typedef DWORD tictoc;
#else
typedef time_t tictoc;
#endif

void tic(tictoc* tv); /* start timing. */

double toc(tictoc* tv); /* stop  timing. */
