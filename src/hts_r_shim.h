#ifndef YANO_HTS_R_SHIM_H
#define YANO_HTS_R_SHIM_H

/*
 * Force-included when compiling the vendored htslib (see Makevars.in).
 * R's compiled-code check rejects packages whose shared object references
 * libc entry points such as abort/exit/sprintf/stdout/stderr.  The macros
 * below redirect those identifiers here, so htslib reports failures
 * through R's error mechanism and writes diagnostics to the process
 * streams opened by yano_shim_init().
 *
 * System headers are included before the macros are installed because
 * macOS <stdio.h> defines its own fortified sprintf/stdout/stderr
 * macros, which would otherwise override definitions supplied on the
 * compiler command line.
 */

#include <stdio.h>
#include <stdlib.h>

#if defined(__GNUC__) || defined(__clang__)
#define YANO_NORET __attribute__((noreturn))
#else
#define YANO_NORET
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern FILE *yano_stdout;
extern FILE *yano_stderr;

int yano_sprintf(char *str, const char *format, ...);
void yano_abort(void) YANO_NORET;
void yano_exit(int status) YANO_NORET;
void yano_shim_init(void);

#ifdef __cplusplus
}
#endif

#undef sprintf
#define sprintf yano_sprintf

#undef stdout
#define stdout yano_stdout

#undef stderr
#define stderr yano_stderr

#undef abort
#define abort yano_abort

#undef exit
#define exit yano_exit

#endif /* YANO_HTS_R_SHIM_H */
