
// config.h

/* Define to dummy 'main' function (if any) required to link to the Fortran libraries. */
//#define FC_DUMMY_MAIN main

/* Define if F77 and FC dummy 'main' functions are identical. */
//#define FC_DUMMY_MAIN_EQ_F77 1

/* Define to a macro mangling the given C identifier (in lower and uppercase), which must not contain underscores, for linking with Fortran. */
//#define FC_FUNC(name, NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
//#define FC_FUNC_(name, NAME) name ## _

/* Define to alternate name for 'main' routine that is called from a 'main' in the Fortran libraries. */
//#define FC_MAIN main

/* Define to 1 if your system has the clock_gettime function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if Fortran interface is to be compiled. */
#define HAVE_FORTRAN 1

/* Define to 1 if you have the <fpu_control.h> header file. */
#define HAVE_FPU_CONTROL_H 1

/* Define to 1 if you have the 'gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <ieeefp.h> header file. */
#define HAVE_IEEEFP_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the 'm' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if stdbool.h conforms to C99. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type '_Bool'. */
#define HAVE__BOOL 1

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
//#define LT_OBJDIR ".libs/"

/* qd major version number */
#define MAJOR_VERSION 2

/* qd minor version number */
#define MINOR_VERSION 3

/* Name of package */
#define PACKAGE "qd"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "yozo@cs.berkeley.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "qd"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "qd 2.3.20"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "qd"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.3.20"

/* qd patch number (sub minor version) */
#define PATCH_VERSION 12

/* Any special symbols needed for exporting APIs. */
//#define QD_API

/* Define this macro to be the copysign(x, y) function. */
#define QD_COPYSIGN(x, y) std::copysign(x, y)

/* Define to 1 to enable debugging code. */
#define QD_DEBUG 1

/* If fused multiply-add is available, define correct macro for using it. */
//#define QD_FMA

/* If fused multiply-subtract is available, define correct macro for using it. */
//#define QD_FMS

/* Define to 1 if your compiler have the C++ standard include files. */
#define QD_HAVE_STD 1

/* Define to 1 to use additions with IEEE-style error bounds. */
#define QD_IEEE_ADD 1

/* Define to 1 to inline commonly used functions. */
//#define QD_INLINE 1

/* Define this macro to be the isfinite(x) function. */
#define QD_ISFINITE(x) std::isfinite(x)

/* Define this macro to be the isinf(x) function. */
#define QD_ISINF(x) std::isinf(x)

/* Define this macro to be the isnan(x) function. */
#define QD_ISNAN(x) std::isnan(x)

/* Define to 1 to use sloppy division (which is faster by slightly inaccurate). */
//#define QD_SLOPPY_DIV 1

/* Define to 1 to use sloppy multiplication (which is faster by slightly inaccurate). */
//#define QD_SLOPPY_MUL 1

/* Set to 1 if using VisualAge C++ compiler for __fmadd builtin. */
//#define QD_VACPP_BUILTINS_H 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define to 1 if your <sys/time.h> declares 'struct tm'. */
#define TM_IN_SYS_TIME 1

/* Version number of package */
#define VERSION "2.3.20"

/* Whether to use x86 fpu fix. */
//#define X86 1
