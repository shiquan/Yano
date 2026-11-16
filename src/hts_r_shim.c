#include <R.h>
#include <R_ext/Error.h>
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>

#include "hts_r_shim.h"

FILE *yano_stdout = NULL;
FILE *yano_stderr = NULL;

int yano_sprintf(char *str, const char *format, ...)
{
    va_list ap;
    int n;

    va_start(ap, format);
    n = vsnprintf(str, (size_t)-1, format, ap);
    va_end(ap);

    return n;
}

void yano_abort(void)
{
    Rf_error("htslib aborted: internal assertion or invariant failure");
}

void yano_exit(int status)
{
    Rf_error("htslib called exit(%d)", status);
}

/* FILE streams for htslib diagnostics.  The process file descriptors are
   wrapped here because R_Consolefile/R_Outputfile are not part of the
   public API. */
static FILE *yano_open_fd(int fd, const char *mode)
{
#ifdef _WIN32
    return _fdopen(fd, mode);
#else
    return fdopen(fd, mode);
#endif
}

void yano_shim_init(void)
{
    if (yano_stdout == NULL)
        yano_stdout = yano_open_fd(1, "w");

    if (yano_stderr == NULL) {
        yano_stderr = yano_open_fd(2, "w");
        if (yano_stderr != NULL)
            setvbuf(yano_stderr, NULL, _IONBF, 0);
    }
}
