#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>

#define ANSI_COLOR_RED     ""
#define ANSI_COLOR_GREEN   ""
#define ANSI_COLOR_YELLOW  ""
#define ANSI_COLOR_BLUE    ""
#define ANSI_COLOR_MAGENTA ""
#define ANSI_COLOR_CYAN    ""
#define ANSI_COLOR_GRAY    ""
#define ANSI_COLOR_RESET   ""

#ifndef error
#define error(line, ...) do						\
    {									\
	Rprintf(ANSI_COLOR_RED "[error] [func: %s, line: %d] " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA line ANSI_COLOR_RESET "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
	errno = 0;							\
	exit(EXIT_FAILURE);						\
    } while(0)

#endif

#define warnings(line, ...) do						\
    {									\
        Rprintf(ANSI_COLOR_YELLOW "[warnings] " line ANSI_COLOR_RESET "\n", ##__VA_ARGS__); \
    } while(0)

#define debug_print(line, ...) do {\
	Rprintf("[ ** DEBUG ** func: %s, line: %d ] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while(0)

#define LOG_print(line, ...) do {               \
	time_t second;                          \
	time(&second);                                                  \
	char _time_buff[100];                                           \
	strftime (_time_buff, 100, "%Y-%m-%d %H:%M:%S", localtime (&second)); \
	Rprintf("[%s] " ANSI_COLOR_GREEN line ANSI_COLOR_RESET"\n", _time_buff, ##__VA_ARGS__); \
    } while(0)

// timer, copied from bwa/utils.c
#include <sys/resource.h>
#include <sys/time.h>

static inline double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
// copied from minimap2/misc.c
static inline long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

#endif
