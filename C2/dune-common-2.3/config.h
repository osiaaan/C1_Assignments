/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.3"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 3

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* does the compiler support __attribute__((deprecated))? */
#define HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#define HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#define HAS_ATTRIBUTE_UNUSED 1

/* Define to 1 if the <array> C++0x is available and support array::fill */
#define HAVE_ARRAY 1

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* define if the Boost library is available */
/* #undef HAVE_BOOST */

/* Define to 1 if you have <boost/make_shared.hpp>. */
/* #undef HAVE_BOOST_MAKE_SHARED_HPP */

/* Define to 1 if you have the <boost/shared_ptr.hpp> header file. */
/* #undef HAVE_BOOST_SHARED_PTR_HPP */

/* Define to 1 if C++11 constexpr is supported */
#define HAVE_CONSTEXPR 1

/* does the compiler support abi::__cxa_demangle */
#define HAVE_CXA_DEMANGLE 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to ENABLE_BOOST if the Boost library is available */
/* #undef HAVE_DUNE_BOOST */

/* Was GMP found and GMP_CPPFLAGS used? */
/* #undef HAVE_GMP */

/* Define to 1 if std::initializer_list is supported */
#define HAVE_INITIALIZER_LIST 1

/* Define to 1 if std::integral_constant< T, v > is supported and casts into T
   */
#define HAVE_INTEGRAL_CONSTANT 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define to 1 if you have the `gmp' library (-L$with_gmp/lib -lgmp). */
/* #undef HAVE_LIBGMP */

/* Define to 1 if you have the `gmpxx' library (-L$with_gmp/lib -lgmpxx). */
/* #undef HAVE_LIBGMPXX */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if SHARED_PTR_NAMESPACE::make_shared is usable. */
#define HAVE_MAKE_SHARED 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <memory> header file. */
#define HAVE_MEMORY 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have the MPI library. This is only true if MPI was found by
   configure _and_ if the application uses the DUNEMPICPPFLAGS (or the
   deprecated MPI_CPPFLAGS) */
/* #undef HAVE_MPI */

/* Define to 1 if you have the symbol mprotect. */
#define HAVE_MPROTECT 1

/* Define to 1 if nullptr is supported */
#define HAVE_NULLPTR 1

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#define HAVE_RPC_RPC_H 1

/* Define to 1 if rvalue references are supported */
#define HAVE_RVALUE_REFERENCES 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if static_assert is supported */
#define HAVE_STATIC_ASSERT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if C++11 std::conditional is supported */
#define HAVE_STD_CONDITIONAL 1

/* Define to 1 if the std::hash template from C++11 is available */
#define HAVE_STD_HASH 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have <sys/mman.h>. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the std::tr1::hash template from TR1 is available */
#define HAVE_TR1_HASH 1

/* Define to 1 if you have the <tr1/memory> header file. */
#define HAVE_TR1_MEMORY 1

/* Define to 1 if you have the <tr1/tuple> header file. */
#define HAVE_TR1_TUPLE 1

/* Define to 1 if you have the <tr1/type_traits> header file. */
#define HAVE_TR1_TYPE_TRAITS 1

/* Define to 1 if you have the <tuple> header file. */
#define HAVE_TUPLE 1

/* Define to 1 if you have the <type_traits> header file. */
#define HAVE_TYPE_TRAITS 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if variadic templates are supported */
#define HAVE_VARIADIC_TEMPLATES 1

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Define to 1 if MPI supports MPI-2 */
/* #undef MPI_2 */

/* Name of package */
#define PACKAGE "dune-common"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "dune-devel@dune-project.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "dune-common"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "dune-common 2.3"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "dune-common"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.3"

/* The header in which SHARED_PTR can be found */
#define SHARED_PTR_HEADER <memory>

/* The namespace in which SHARED_PTR can be found */
#define SHARED_PTR_NAMESPACE std

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Version number of package */
#define VERSION "2.3"

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

#include <dune/common/deprecated.hh>

#include <dune/common/unused.hh>
