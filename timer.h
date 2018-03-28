/* This unit contains three functions for measuring cpu time	*/
/* used by a program.  The first, user_time(), returns the	*/
/* total amount of user cpu time used by the current process.	*/
/* The second, start_timer(), starts the timer.			*/
/* Actually, it only records the current value of user_time()	*/
/* for future reference.  Finally, elapsed_time() returns the	*/
/* amount of cpu time that has been used since the last call	*/
/* to start_timer().						*/

	/* ERROR_VALUE is returned by elapsed_time()	*/
	/* if start_time() was not called previously.	*/
#define  ERROR_VALUE -1.0

	/* Return the amount of time in seconds used by	*/
	/* the current process since it began.		*/
#ifdef __STDC__
double user_time(void);
#else
double user_time();
#endif


	/* Starts timer. */
#ifdef __STDC__
void start_timer(void);
#else
start_timer();
#endif


	/* Returns elapsed time since last call to start_timer() */
	/* Returns ERROR_VALUE if start_timer() has never been called */
#ifdef __STDC__
double elapsed_time(void);
#else
double elapsed_time();
#endif
