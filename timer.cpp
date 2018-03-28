#include <sys/time.h>
#include <sys/resource.h>
#include "timer.h"

#define FALSE 0
#define TRUE 1


static int  timer_set = FALSE;
static double   old_time;


/* Return the amount of time in seconds used by the current process since
   it began. */
#ifdef __STDC__
double  user_time (void)
#else
double  user_time ()
#endif
{
    struct rusage   t;

/* Get resource usage from system */
    getrusage (RUSAGE_SELF, &t);
/* Return processes userr mode time in seconds. */
    return (t.ru_utime.tv_sec + (t.ru_utime.tv_usec / 1000000.0));
}


/* Starts timer. */
#ifdef __STDC__
void start_timer (void)
#else
start_timer ()
#endif
{
    timer_set = TRUE;
    old_time = user_time ();
}


/* Returns elapsed time since last call to start_timer().  Returns ERROR_VALUE */
/* if Start_Timer() has never been called. */
#ifdef __STDC__
double  elapsed_time (void)
#else
double  elapsed_time ()
#endif
{
    if      (timer_set == FALSE) {
	        return (ERROR_VALUE);
    }
    else {
	return (user_time () - old_time);
    }
}


#ifdef DEBUG
#define MAX_ITERATIONS 100000

main (argc, argv)
int     argc;
char   *argv[];
{
    int     i = argc == 1 ? MAX_ITERATIONS : atoi (argv[1]);

    start_timer ();
    while (--i);
    printf ("time = %f\n", elapsed_time ());
}
#endif
