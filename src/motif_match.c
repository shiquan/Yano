#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"
#include "_XVector_stubs.c"

static ByteTrTable byte2offset;

static int no_warning_yet;

static void set_byte2offset_elt(ByteTrTable *byte2offset, int byte, int offset, int error_on_dup)
{
	int *offset_p;

	if (byte < 0 || byte >= BYTETRTABLE_LENGTH)
		error("Biostrings internal error in set_byte2offset_elt(): "
		      "invalid byte value %d", byte);
	offset_p = byte2offset->byte2code + (unsigned char) byte;
	if (*offset_p == NA_INTEGER) {
		*offset_p = offset;
		return;
	}
	if (error_on_dup)
		error("Biostrings internal error in set_byte2offset_elt(): "
		      "duplicated byte value %d", byte);
	return;
}

/*
 * Values in 'bytes' must represent byte values i.e. values >= 0 and < 256.
 * The byte offsets are written to 'byte2offset'.
*/
void _init_byte2offset_with_INTEGER(ByteTrTable *byte2offset, SEXP bytes, int error_on_dup)
{
    int byte, offset;
    
    if (LENGTH(bytes) > BYTETRTABLE_LENGTH)
        error("Biostrings internal error in "
              "_init_byte2offset_with_INTEGER(): "
              "LENGTH(bytes) > BYTETRTABLE_LENGTH");
    for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
        byte2offset->byte2code[byte] = NA_INTEGER;
    for (offset = 0; offset < LENGTH(bytes); offset++) {
        byte = INTEGER(bytes)[offset];
        set_byte2offset_elt(byte2offset, byte, offset, error_on_dup);
    }
    return;
}

static double compute_pwm_score(const double *pwm, int pwm_ncol, const char *S, int nS, int pwm_shift)
{
    int i, rowoffset;
    double score;
    
    S += pwm_shift;
    nS -= pwm_shift;
    if (pwm_shift < 0 || nS < pwm_ncol)
        error("'starting.at' contains invalid values");
    score = 0.00;
    for (i = 0; i < pwm_ncol; i++, pwm += 4, S++) {
        rowoffset = byte2offset.byte2code[(unsigned char) *S];
        if (rowoffset == NA_INTEGER) {
            if (no_warning_yet) {
                warning("'subject' contains letters not in "
                        "[ACGT] ==> assigned weight 0 to them");
                no_warning_yet = 0;
            }
            continue;
        }
        score += pwm[rowoffset];
    }
    return score;
}

static int _match_PWM(const double *pwm, int pwm_ncol, const Chars_holder *S, double minscore)
{
    int n1, n2;
    double score;
    for (n1 = 0, n2 = pwm_ncol; n2 <= S->length; n1++, n2++) {
        score = compute_pwm_score(pwm, pwm_ncol, S->ptr, S->length, n1);
        if (score >= minscore) return 1;
    }
    return 0;
}

SEXP match_PWM_fast(SEXP pwm, SEXP subject, SEXP min_score, SEXP count_only, SEXP base_codes)
{
    Chars_holder S;
    int pwm_ncol, is_count_only;
    double minscore;
    
    if (INTEGER(GET_DIM(pwm))[0] != 4)
        error("'pwm' must have 4 rows");
    pwm_ncol = INTEGER(GET_DIM(pwm))[1];
    S = hold_XRaw(subject);

    Rprintf("%d\t%s\n", S.length, S.ptr);
    minscore = REAL(min_score)[0];
    
    _init_byte2offset_with_INTEGER(&byte2offset, base_codes, 1);
    no_warning_yet = 1;

    int ret = _match_PWM(REAL(pwm), pwm_ncol, &S, minscore);

    return ScalarInteger(ret);
}

    
