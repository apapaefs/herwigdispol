/* Config/config.h.  Generated from config.h.in by configure.  */
/* Config/config.h.in.  Generated from configure.ac by autoheader.  */

/* Defined if the requested minimum BOOST version is satisfied */
#define HAVE_BOOST 1

/* Define to 1 if you have <boost/test/unit_test.hpp> */
#define HAVE_BOOST_TEST_UNIT_TEST_HPP 1

/* Defined if the Boost unit_test_framework library is available */
#define HAVE_BOOST_UNIT_TEST_FRAMEWORK 1

/* define if the compiler supports basic C++11 syntax */
#define HAVE_CXX11 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* We have HepMC3 */
#define HAVE_HEPMC3 1

/* Define to 1 if you have the <HepMC3/WriterRootTree.h> header file. */
/* #undef HAVE_HEPMC3_WRITERROOTTREE_H */

/* Define to 1 if you have the <HepMC3/WriterRoot.h> header file. */
/* #undef HAVE_HEPMC3_WRITERROOT_H */

/* Define to 1 if you have the <HepMC3/Writer.h> header file. */
#define HAVE_HEPMC3_WRITER_H 1

/* Define to 1 if you have the <HepMC/HepMCDefs.h> header file. */
/* #undef HAVE_HEPMC_HEPMCDEFS_H */

/* Define to 1 if you have the <HepMC/PdfInfo.h> header file. */
/* #undef HAVE_HEPMC_PDFINFO_H */

/* Define to 1 if you have the <history.h> header file. */
/* #undef HAVE_HISTORY_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `gsl' library (-lgsl). */
/* #undef HAVE_LIBGSL */

/* Define to 1 if you have the `gslcblas' library (-lgslcblas). */
/* #undef HAVE_LIBGSLCBLAS */

/* Define to 1 if you have the `m' library (-lm). */
/* #undef HAVE_LIBM */

/* Define if you have a readline compatible library */
/* #undef HAVE_LIBREADLINE */

/* Define to 1 if you have the `Rivet' library (-lRivet). */
/* #undef HAVE_LIBRIVET */

/* Define to 1 if you have `z' library (-lz) */
#define HAVE_LIBZ 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <readline.h> header file. */
/* #undef HAVE_READLINE_H */

/* Define if your readline library has \`add_history' */
/* #undef HAVE_READLINE_HISTORY */

/* Define to 1 if you have the <readline/history.h> header file. */
/* #undef HAVE_READLINE_HISTORY_H */

/* Define to 1 if you have the <readline/readline.h> header file. */
/* #undef HAVE_READLINE_READLINE_H */

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

/* Has GenCrossection */
#define HEPMC_HAS_CROSS_SECTION 1

/* Has named weights */
#define HEPMC_HAS_NAMED_WEIGHTS 1

/* Has GenPdfInfo */
#define HEPMC_HAS_PDF_INFO 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "ThePEG"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://www.thep.lu.se/ThePEG/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ThePEG"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ThePEG devel"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "ThePEG"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "devel"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* The command which, taking the name of a bzipped file as argument, unzips it
   and prints it to stdout. Default is "bunzip2 -c". */
#define ThePEG_BZ2READ_FILE "bunzip2 -c"

/* The command which, taking the name of a bzipped file as argument, reads
   stdin, zips it and writes it to the file. Default is "bzip2 -c > ". */
#define ThePEG_BZ2WRITE_FILE "bzip2 -c > "

/* The command which, taking the name of a gzipped file as argument, unzips it
   and prints it to stdout. Default is "gunzip -c" */
#define ThePEG_GZREAD_FILE "gunzip -c"

/* The command which, taking the name of a gzipped file as argument, reads
   stdin, zips it and writes it to the file. Default is "gzip -c > ". */
#define ThePEG_GZWRITE_FILE "gzip -c > "

/* define if dlopen is available */
#define ThePEG_HAS_DLOPEN 1

/* define if expm1 is available */
#define ThePEG_HAS_EXPM1 1

/* define if fenv is available */
#define ThePEG_HAS_FENV 1

/* define if fpucontrol is available */
#define ThePEG_HAS_FPU_CONTROL 1

/* define if log1p is available */
#define ThePEG_HAS_LOG1P 1

/* Rivet major version (1,2,3) */
#define ThePEG_RIVET_VERSION 4

/* Version number of package */
#define VERSION "devel"
