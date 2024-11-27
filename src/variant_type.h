#ifndef SEQUENCE_HEADER
#define SEQUENCE_HEADER

#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER) && !defined(__clang__)
# define inline __inline
#endif

#define C4_A  0
#define C4_C  1
#define C4_G  2
#define C4_T  3
#define C4_U  3
#define C4_N  4
#define seqarr    "ACGTN"
#define revseqarr "TGCAN"

#define SEQ_COMP(a,b) (a + b == 3)

typedef char * (*func_dup_seq)(const char *, unsigned long );

extern char *rev_seqs(const char *dna_seqs, unsigned long n);
extern int check_stop_codon(char *seq, char *p_end, int mito);
extern void compl_seq(char *seq, int l);
extern int seq2code4(int seq);
extern int same_DNA_seqs(const char *a, const char *b, int l );
extern int codon2aminoid(char *codon,int mito);
extern char *rev_seqs(const char *dna_seqs, unsigned long n);

#define X_CODO   0

#define C4_Stop  0
#define C4_Phe   1
#define C4_Leu   2
#define C4_Ser   3
#define C4_Tyr   4
#define C4_Cys   5
#define C4_Trp   6
#define C4_Pro   7
#define C4_His   8
#define C4_Gln   9
#define C4_Arg  10
#define C4_Ile  11
#define C4_Met  12
#define C4_Thr  13
#define C4_Asn  14
#define C4_Lys  15
#define C4_Val  16
#define C4_Ala  17
#define C4_Asp  18
#define C4_Glu  19
#define C4_Gly  20
#define C4_Sec  21 // Selenocysteine

static inline char *codon2name(int code) {
    static char *codon_names[] = {
        "*", "Phe", "Leu", "Ser", "Tyr", "Cys", "Trp", "Pro", "His", "Gln", "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Val", "Ala", "Asp", "Glu", "Gly", "Sec"
};
    return codon_names[code];
}
static inline char *codon2name2(int code) {
    static char *codon_short_names[] = { "*", "F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G", "U"};
    return codon_short_names[code];
}
const static int codon_matrix[4][4][4] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
      { C4_Ile, C4_Ile, C4_Met, C4_Ile, }, },
    
    { { C4_Gln, C4_His, C4_Gln, C4_His, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, }, },
    
    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
      { C4_Val, C4_Val, C4_Val, C4_Val, }, },
    
    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
      { C4_Stop, C4_Cys, C4_Trp, C4_Cys, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, }, },
};

/*
  https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoCode

  For human mitochondrial genes, unlike the universal code, UGA codes for tryptophan instead of termination and AUA codes for methionine instead of isoleucine.
 */

const static int mitomap_codon_matrix[4][4][4] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
      { C4_Met, C4_Ile, C4_Met, C4_Ile, }, },
    
    { { C4_Gln, C4_His, C4_Gln, C4_His, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, }, },
    
    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
      { C4_Val, C4_Val, C4_Val, C4_Val, }, },
    
    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
      { C4_Trp, C4_Cys, C4_Trp, C4_Cys, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, }, },
};

// 1 on yes, 0 on no
static inline int check_is_stop(char *codon, int mito)
{
    if ( mito == 0 ) 
        return codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
    else
        return mitomap_codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
}

#endif
