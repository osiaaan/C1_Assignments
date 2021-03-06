/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to the version of C2-practical */
#define C2_PRACTICAL_VERSION "1"

/* Define to the major version of C2-practical */
#define C2_PRACTICAL_VERSION_MAJOR 1

/* Define to the minor version of C2-practical */
#define C2_PRACTICAL_VERSION_MINOR 0

/* Define to the revision of C2-practical */
#define C2_PRACTICAL_VERSION_REVISION 0

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0
   */
/* #undef DUNE_ALBERTA_VERSION */

/* Define to the version of dune-alugrid */
#define DUNE_ALUGRID_VERSION "2.3"

/* Define to the major version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MAJOR 2

/* Define to the minor version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MINOR 3

/* Define to the revision of dune-alugrid */
#define DUNE_ALUGRID_VERSION_REVISION 0

/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.3"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 3

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "2.3"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR 3

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION 0

/* If this is set, public access to the implementation of facades like Entity,
   Geometry, etc. is granted. */
/* #undef DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS */

/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "2.3"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR 2

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR 3

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION 0

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

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
/* #undef HAVE_ALBERTA */

/* This is only true if alugrid-library was found by configure _and_ if the
   application uses the ALUGRID_CPPFLAGS */
/* #undef HAVE_ALUGRID */

/* Define to 1 if you have the <alugrid_parallel.h> header file. */
/* #undef HAVE_ALUGRID_PARALLEL_H */

/* Define to 1 if you have the <alugrid_serial.h> header file. */
/* #undef HAVE_ALUGRID_SERIAL_H */

/* Define to 1 if amiramesh-library is found */
/* #undef HAVE_AMIRAMESH */

/* Use the Apple OpenGL framework. */
/* #undef HAVE_APPLE_OPENGL_FRAMEWORK */

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

/* Define to 1 if dune-alugrid was found */
#define HAVE_DUNE_ALUGRID 1

/* Define to ENABLE_BOOST if the Boost library is available */
/* #undef HAVE_DUNE_BOOST */

/* Define to 1 if dune-common was found */
#define HAVE_DUNE_COMMON 1

/* Define to 1 if dune-geometry was found */
#define HAVE_DUNE_GEOMETRY 1

/* Define to 1 if dune-grid was found */
#define HAVE_DUNE_GRID 1

/* Was GMP found and GMP_CPPFLAGS used? */
/* #undef HAVE_GMP */

/* This is only true if grape-library was found by configure _and_ if the
   application uses the GRAPE_CPPFLAGS */
/* #undef HAVE_GRAPE */

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

/* Define if you have METIS library */
/* #undef HAVE_METIS */

/* Define if you have the MPI library. This is only true if MPI was found by
   configure _and_ if the application uses the DUNEMPICPPFLAGS (or the
   deprecated MPI_CPPFLAGS) */
/* #undef HAVE_MPI */

/* Define to 1 if you have the symbol mprotect. */
#define HAVE_MPROTECT 1

/* Define to 1 if nullptr is supported */
#define HAVE_NULLPTR 1

/* Define if you have the Parmetis library. This is only true if MPI was found
   by configure _and_ if the application uses the PARMETIS_CPPFLAGS */
/* #undef HAVE_PARMETIS */

/* Define to 1 if psurface-library is found */
/* #undef HAVE_PSURFACE */

/* If set we have at least psurface version 2.0 */
/* #undef HAVE_PSURFACE_2_0 */

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

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

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
/* #undef HAVE_UG */

/* Do we have UG in at least version 3.9.1-patch10? */
/* #undef HAVE_UG_PATCH10 */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if variadic templates are supported */
#define HAVE_VARIADIC_TEMPLATES 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* This is only true if zoltan-library was found by configure _and_ if the
   application uses the ZOLTAN_CPPFLAGS */
/* #undef HAVE_ZOLTAN */

/* Define to 1 if you have the <zoltan_cpp.h> header file. */
/* #undef HAVE_ZOLTAN_CPP_H */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Define to 1 if MPI supports MPI-2 */
/* #undef MPI_2 */

/* Name of package */
#define PACKAGE "c2-practical"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "A.S.Dedner@warwick.ac.uk"

/* Define to the full name of this package. */
#define PACKAGE_NAME "C2-practical"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "C2-practical 1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "c2-practical"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1"

/* The namespace prefix of the psurface library */
/* #undef PSURFACE_NAMESPACE */

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* The header in which SHARED_PTR can be found */
#define SHARED_PTR_HEADER <memory>

/* The namespace in which SHARED_PTR can be found */
#define SHARED_PTR_NAMESPACE std

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Version number of package */
#define VERSION "1"

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */

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

/* add GRIDTYPE typedef for grid implementation Dune::YaspGrid< dimgrid >:
    defining YASPGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined YASPGRID && ! defined USED_YASPGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (GRIDDIM == WORLDDIM)
    #error "Preprocessor assertion GRIDDIM == WORLDDIM failed."
  #endif

  #include <dune/grid/yaspgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfyasp.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::YaspGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_YASPGRID_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined YASPGRID && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::conforming >:
    defining ALUGRID_CONFORM during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_CONFORM && ! defined USED_ALUGRID_CONFORM_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/alugrid/grid.hh>
  #include <dune/alugrid/dgf.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::conforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CONFORM_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_CONFORM && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, Dune::cube, Dune::nonconforming >:
    defining ALUGRID_CUBE during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_CUBE && ! defined USED_ALUGRID_CUBE_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/alugrid/grid.hh>
  #include <dune/alugrid/dgf.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, Dune::cube, Dune::nonconforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CUBE_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_CUBE && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::nonconforming >:
    defining ALUGRID_SIMPLEX during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_SIMPLEX && ! defined USED_ALUGRID_SIMPLEX_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/alugrid/grid.hh>
  #include <dune/alugrid/dgf.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::nonconforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_SIMPLEX_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_SIMPLEX && ..

#include <dune/common/unused.hh>

/* add GRIDTYPE typedef for grid implementation Dune::AlbertaGrid< dimgrid >:
    defining ALBERTAGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALBERTAGRID && ! defined USED_ALBERTAGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (WORLDDIM == ALBERTA_DIM)
    #error "Preprocessor assertion WORLDDIM == ALBERTA_DIM failed."
  #endif

  #include <dune/grid/albertagrid.hh>
  #include <dune/grid/albertagrid/dgfparser.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::AlbertaGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALBERTAGRID_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALBERTAGRID && ..

/* add GRIDTYPE typedef for grid implementation Dune::UGGrid< dimgrid >:
    defining UGGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined UGGRID && ! defined USED_UGGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! (GRIDDIM == WORLDDIM)
    #error "Preprocessor assertion GRIDDIM == WORLDDIM failed."
  #endif

  #include <dune/grid/uggrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfug.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::UGGrid< dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_UGGRID_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined UGGRID && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, simplex, conforming >:
    defining ALUGRID_CONFORM during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_CONFORM && ! defined USED_ALUGRID_CONFORM_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, simplex, conforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CONFORM_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_CONFORM && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, cube, nonconforming >:
    defining ALUGRID_CUBE during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_CUBE && ! defined USED_ALUGRID_CUBE_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, cube, nonconforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_CUBE_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_CUBE && ..

/* add GRIDTYPE typedef for grid implementation Dune::ALUGrid< dimgrid, dimworld, simplex, nonconforming >:
    defining ALUGRID_SIMPLEX during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ALUGRID_SIMPLEX && ! defined USED_ALUGRID_SIMPLEX_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/alugrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfalu.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::ALUGrid< dimgrid, dimworld, simplex, nonconforming > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ALUGRID_SIMPLEX_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ALUGRID_SIMPLEX && ..

/* add GRIDTYPE typedef for grid implementation Dune::OneDGrid:
    defining ONEDGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined ONEDGRID && ! defined USED_ONEDGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #if ! ((GRIDDIM == 1) && (WORLDDIM == 1))
    #error "Preprocessor assertion (GRIDDIM == 1) && (WORLDDIM == 1) failed."
  #endif

  #include <dune/grid/onedgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfoned.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::OneDGrid GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_ONEDGRID_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined ONEDGRID && ..

/* add GRIDTYPE typedef for grid implementation Dune::SGrid< dimgrid, dimworld >:
    defining SGRID during compilation typedefs this grid implementation as GridType
    in namespace Dune::GridSelector;
    also integer constants dimgrid and dimworld are set in this namespace.
    The required headers for this grid implementation are also included.
  */
 #if HAVE_DUNE_GRID && defined SGRID && ! defined USED_SGRID_GRIDTYPE
  #if HAVE_GRIDTYPE
   #error "Ambiguous definition of GRIDTYPE."
  #endif 

  #ifndef WORLDDIM
    #define WORLDDIM GRIDDIM
  #endif
  #if not (WORLDDIM >= GRIDDIM)
    #error "WORLDDIM < GRIDDIM does not make sense."
  #endif

  #include <dune/grid/sgrid.hh>
  #include <dune/grid/io/file/dgfparser/dgfs.hh>

  namespace Dune
  {
    namespace GridSelector
    {
      const int dimgrid = GRIDDIM;
      const int dimworld = WORLDDIM;
      typedef Dune::SGrid< dimgrid, dimworld > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
  #define USED_SGRID_GRIDTYPE 1
#endif // #if HAVE_DUNE_GRID && defined SGRID && ..
