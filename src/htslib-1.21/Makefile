# Makefile for htslib — minimal static build for bundling in R packages
#
# Copyright (C) 2013-2024 Genome Research Ltd.
# Adapted for minimal dependency builds.
#
# Usage:
#   make lib-static                          # build libhts.a with default settings
#   make lib-static USE_BZ2=0 USE_LZMA=0 USE_CURL=0  # minimal deps (zlib only)

CC     ?= gcc
AR     ?= ar
RANLIB ?= ranlib

# ---- optional library auto-detection ----
# Set to 0 on the command line to disable, or 1 to force-enable.
# Auto-detect probes try to compile+link a trivial program.
USE_BZ2 ?= $(shell echo 'int main(){return 0;}' | $(CC) -include bzlib.h -x c - -lbz2 -o /dev/null 2>/dev/null && echo 1 || echo 0)
USE_LZMA ?= $(shell echo 'int main(){return 0;}' | $(CC) -include lzma.h -x c - -llzma -o /dev/null 2>/dev/null && echo 1 || echo 0)
USE_CURL ?= $(shell echo 'int main(){return 0;}' | $(CC) -include curl/curl.h -x c - -lcurl -o /dev/null 2>/dev/null && echo 1 || echo 0)

# ---- required flags ----
CFLAGS   ?= -g -O2 -fpic
CPPFLAGS := -I. $(CPPFLAGS)
LDFLAGS  ?= -fpic

# ---- platform detection ----
PLATFORM := $(shell uname -s)

# ---- htslib configuration ----
htslib_default_libs = -lz -lm
ifeq ($(USE_BZ2),1)
  CPPFLAGS += -DHAVE_LIBBZ2
  htslib_default_libs += -lbz2
endif
ifeq ($(USE_LZMA),1)
  CPPFLAGS += -DHAVE_LIBLZMA
  ifneq ($(PLATFORM),Darwin)
    CPPFLAGS += -DHAVE_LZMA_H
  endif
  htslib_default_libs += -llzma
endif
ifeq ($(USE_CURL),1)
  CPPFLAGS += -DHAVE_LIBCURL
  htslib_default_libs += -lcurl
  NONCONFIGURE_OBJS = hfile_libcurl.o
else
  NONCONFIGURE_OBJS =
endif

# ---- SIMD: auto-detect via actual source compilation ----
# Probe by compiling the real htscodecs SIMD source with the ISA flag.
# If the bundled htscodecs is incompatible with the system compiler
# (e.g. GCC 13 ABI changes), the extension is silently disabled.
ARCH := $(shell uname -m)

HTS_CFLAGS_SSE4  = -msse4.1
HTS_CFLAGS_AVX2  = -mavx2
HTS_CFLAGS_AVX512 = -mavx512f -mavx512bw

ifneq ($(filter x86_64 amd64,$(ARCH)),)
  HTS_BUILD_SSE4 := $(shell $(CC) $(HTS_CFLAGS_SSE4) -I. -c $(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_sse4.c -o /dev/null 2>/dev/null && echo 1 || echo 0)
  HTS_BUILD_AVX2 := $(shell $(CC) $(HTS_CFLAGS_AVX2) -I. -c $(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_avx2.c -o /dev/null 2>/dev/null && echo 1 || echo 0)
  HTS_BUILD_AVX512 := $(shell $(CC) $(HTS_CFLAGS_AVX512) -I. -c $(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_avx512.c -o /dev/null 2>/dev/null && echo 1 || echo 0)
  ifeq ($(HTS_BUILD_SSE4),0)
    HTS_CFLAGS_SSE4 =
  endif
  ifeq ($(HTS_BUILD_AVX2),0)
    HTS_CFLAGS_AVX2 =
  endif
  ifeq ($(HTS_BUILD_AVX512),0)
    HTS_CFLAGS_AVX512 =
  endif
else
  HTS_CFLAGS_SSE4  =
  HTS_CFLAGS_AVX2  =
  HTS_CFLAGS_AVX512 =
  HTS_BUILD_SSE4  =
  HTS_BUILD_AVX2  =
  HTS_BUILD_AVX512 =
endif

# Always enable CPUID check for x86 (compile-guarded by preprocessor)
ifeq ($(filter x86_64 amd64,$(ARCH)),)
  HTS_HAVE_CPUID =
else
  HTS_HAVE_CPUID = 1
endif

LIBS = $(htslib_default_libs)

srcdir = .
srcprefix =
HTSPREFIX =

# ---- htscodecs (bundled) ----
include htscodecs_bundled.mk

# ---- object files ----
LIBHTS_OBJS = \
	kfunc.o \
	kstring.o \
	bcf_sr_sort.o \
	bgzf.o \
	errmod.o \
	faidx.o \
	header.o \
	hfile.o \
	hts.o \
	hts_expr.o \
	hts_os.o \
	md5.o \
	multipart.o \
	probaln.o \
	realn.o \
	regidx.o \
	region.o \
	sam.o \
	sam_mods.o \
	simd.o \
	synced_bcf_reader.o \
	vcf_sweep.o \
	tbx.o \
	textutils.o \
	thread_pool.o \
	vcf.o \
	vcfutils.o \
	cram/cram_codecs.o \
	cram/cram_decode.o \
	cram/cram_encode.o \
	cram/cram_external.o \
	cram/cram_index.o \
	cram/cram_io.o \
	cram/cram_stats.o \
	cram/mFILE.o \
	cram/open_trace_file.o \
	cram/pooled_alloc.o \
	cram/string_alloc.o \
	$(HTSCODECS_OBJS) \
	$(NONCONFIGURE_OBJS)

# ---- required header dependencies ----
header_h = header.h cram/string_alloc.h cram/pooled_alloc.h $(htslib_khash_h) $(htslib_kstring_h) $(htslib_sam_h)
hfile_internal_h = hfile_internal.h $(htslib_hts_defs_h) $(htslib_hfile_h) $(textutils_internal_h)
hts_internal_h = hts_internal.h $(htslib_hts_h) $(textutils_internal_h)
textutils_internal_h = textutils_internal.h $(htslib_kstring_h)
cram_h = cram/cram.h $(cram_samtools_h) $(header_h) $(cram_structs_h) $(cram_io_h) cram/cram_encode.h cram/cram_decode.h cram/cram_stats.h cram/cram_codecs.h cram/cram_index.h $(htslib_cram_h)
cram_io_h = cram/cram_io.h $(cram_misc_h)
cram_misc_h = cram/misc.h
cram_os_h = cram/os.h $(htslib_hts_endian_h)
cram_samtools_h = cram/cram_samtools.h $(htslib_sam_h)
cram_structs_h = cram/cram_structs.h $(htslib_thread_pool_h) $(htslib_cram_h) cram/string_alloc.h cram/mFILE.h $(htslib_khash_h)
cram_open_trace_file_h = cram/open_trace_file.h cram/mFILE.h
bcf_sr_sort_h = bcf_sr_sort.h $(htslib_synced_bcf_reader_h) $(htslib_kbitset_h)
fuzz_settings_h = fuzz_settings.h
hts_time_funcs_h = hts_time_funcs.h
sam_internal_h = sam_internal.h $(htslib_sam_h)
thread_pool_internal_h = thread_pool_internal.h $(htslib_thread_pool_h)

# ---- build targets ----
ALL_CPPFLAGS = -I. $(CPPFLAGS)

# ---- per-file CFLAGS for SIMD objects (needed for -mavx2 etc.) ----
htscodecs/htscodecs/rANS_static32x16pr_sse4.o:   CFLAGS += $(HTS_CFLAGS_SSE4)
htscodecs/htscodecs/rANS_static32x16pr_avx2.o:   CFLAGS += $(HTS_CFLAGS_AVX2)
htscodecs/htscodecs/rANS_static32x16pr_avx512.o: CFLAGS += $(HTS_CFLAGS_AVX512)

all: lib-static

.c.o:
	$(CC) $(CFLAGS) $(ALL_CPPFLAGS) -c -o $@ $<

lib-static: libhts.a

version.h:
	echo '#define HTS_VERSION_TEXT "1.21"' > $@

config_vars.h:
	printf '#define HTS_CC "%s"\n#define HTS_CPPFLAGS "%s"\n#define HTS_CFLAGS "%s"\n#define HTS_LDFLAGS "%s"\n#define HTS_LIBS "%s"\n' \
	'$(CC)' \
	'$(CPPFLAGS)' \
	'$(CFLAGS)' \
	'$(LDFLAGS)' \
	'$(LIBS)' > $@

libhts.a: config.h version.h config_vars.h htscodecs.mk $(LIBHTS_OBJS)
	@-rm -f $@
	$(AR) -rc $@ $(LIBHTS_OBJS)
	-$(RANLIB) $@

# ---- generated config ----
# htscodecs.mk: sets up SIMD codec object list and flags
htscodecs.mk:
	@echo '# Default htscodecs.mk generated by Makefile' > $@
	@echo 'include $$(HTSPREFIX)htscodecs_bundled.mk' >> $@
	@echo '' >> $@
	@echo '# Probe results' >> $@
	@echo "HTS_HAVE_CPUID = $(HTS_HAVE_CPUID)" >> $@
	@echo "HTS_BUILD_SSE4 = $(HTS_BUILD_SSE4)" >> $@
	@echo "HTS_BUILD_AVX2 = $(HTS_BUILD_AVX2)" >> $@
	@echo "HTS_BUILD_AVX512 = $(HTS_BUILD_AVX512)" >> $@

config.h:
	@echo '/* config.h generated by Makefile */' > $@
	@echo '#ifndef _XOPEN_SOURCE' >> $@
	@echo '#define _XOPEN_SOURCE 600' >> $@
	@echo '#endif' >> $@
	@echo '#define HAVE_DRAND48 1' >> $@
	@if [ "$(USE_BZ2)" = "1" ]; then \
		echo '#define HAVE_LIBBZ2 1' >> $@ ; \
	fi
	@if [ "$(USE_LZMA)" = "1" ]; then \
		echo '#define HAVE_LIBLZMA 1' >> $@ ; \
		if [ "$(PLATFORM)" != "Darwin" ]; then \
			echo '#define HAVE_LZMA_H 1' >> $@ ; \
		fi \
	fi
	@if [ "$(USE_CURL)" = "1" ]; then \
		echo '#define HAVE_LIBCURL 1' >> $@ ; \
	fi
	@if [ "$(HTS_HAVE_CPUID)" != "" ]; then \
		echo '#define HAVE_DECL___CPUID_COUNT 1' >> $@ ; \
		echo '#define HAVE_DECL___GET_CPUID_MAX 1' >> $@ ; \
	fi
	@if [ "$(HTS_BUILD_SSE4)" != "" ]; then \
		echo '#define HAVE_POPCNT 1' >> $@ ; \
		echo '#define HAVE_SSE4_1 1' >> $@ ; \
		echo '#define HAVE_SSSE3 1' >> $@ ; \
	fi
	@if [ "$(HTS_BUILD_AVX2)" != "" ]; then \
		echo '#define HAVE_AVX2 1' >> $@ ; \
	fi
	@if [ "$(HTS_BUILD_AVX512)" != "" ]; then \
		echo '#define HAVE_AVX512 1' >> $@ ; \
	fi
	@echo '#if defined __x86_64__ || defined __arm__ || defined __aarch64__' >> $@
	@echo '#define HAVE_ATTRIBUTE_CONSTRUCTOR 1' >> $@
	@echo '#endif' >> $@
	@echo '#if (defined(__x86_64__) || defined(_M_X64))' >> $@
	@echo '#define HAVE_ATTRIBUTE_TARGET 1' >> $@
	@echo '#define HAVE_BUILTIN_CPU_SUPPORT_SSSE3 1' >> $@
	@echo '#endif' >> $@

# ---- htslib_vars.mk from upstream (defines HTSCODECS_OBJS, etc.) ----
include htslib_vars.mk

# ---- conditional inclusion of config.mk (if exists from configure) ----
-include config.mk

# ---- clean ----
clean:
	rm -f *.o *.pico cram/*.o cram/*.pico test/*.o test/*.dSYM
	rm -f htscodecs/htscodecs/*.o htscodecs/htscodecs/*.pico
	rm -f htscodecs/tests/*.o
	rm -f hts-object-files
	rm -f libhts.a libhts.so libhts.dylib libhts.*.dylib
	rm -f annot-tsv bgzip htsfile tabix
	rm -f test/hts_endian test/fieldarith test/hfile test/pileup test/sam
	rm -f test/test_bgzf test/test_expr test/test_faidx test/test_kfunc
	rm -f test/test_khash test/test_kstring test/test_mod test/test_view
	rm -f test/test_index test/test-vcf-api test/test-vcf-sweep test/test-bcf-sr
	rm -f test/test_realn test/test-regidx test/test_str2int test/test_time_funcs
	rm -f test/test-bcf-translate test/test-parse-reg test/test_introspection
	rm -f test/test-bcf_set_variant_type test/thrash_threads*
	rm -f test/fuzz/hts_open_fuzzer.o test/pileup_mod test/plugins-dlhts
	rm -f htscodecs/tests/rans4x8 htscodecs/tests/rans4x16pr
	rm -f htscodecs/tests/arith_dynamic htscodecs/tests/tokenise_name3
	rm -f htscodecs/tests/fqzcomp_qual htscodecs/tests/varint
	rm -f test/*.tmp test/*.tmp.* test/faidx/*.tmp*
	rm -f test/longrefs/*.tmp.* test/tabix/*.tmp.* test/bgzf_boundaries/*.tmp.*
	rm -f test/*/FAIL* header-exports.txt shlib-exports-dylib.txt
	rm -f config.h config_vars.h htscodecs.mk version.h htslib_static.mk
	rm -f htslib.pc.tmp *.pc

.PHONY: all lib-static clean
