#
# After changing this file run `autoconf` to produce the 'configure' script.
#
# MPI discovery options for the build process:
#
# - Run with automatic detection: `R CMD INSTALL {pkg}`
#
# - Run with a prefix folder.  Subfolders will be searched for binarys.
#   `R CMD INSTALL {pkg} --configure-args=--with-mpi=/path/to/location`
#
# - Fully specify paths and type.
#   `R CMD INSTALL {pkg} --configure-args="        \\
#    --with-mpi-include=/path/to/mpi_include_dir  \\
#    --with-mpi-libpath=/path/to/mpi_lib_dir      \\
#    --with-mpi-type={OPENMPI,MPICH,MPICH2,LAM}"`

AC_PREREQ([2.69])
AC_INIT([RSiena], 1.1-289)
AC_CONFIG_SRCDIR([src])

# Find R home and set CC/CFLAGS.  This is required for the custom openmp test.
R_HOME=`R RHOME`
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`
LIBS="${PKG_LIBS}"


AC_PROG_CXX
AC_PROG_CC
AC_CHECK_HEADERS([sys/time.h unistd.h])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime gettimeofday memmove memset pow select sqrt])

##
## OpenMP
##
## Based on the autoconf.ac from the `randomForestSRC` package.
##
AC_OPENMP

ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([whether OpenMP will work in a package])
  AC_LANG_CONFTEST([AC_LANG_PROGRAM(
    [[#include <omp.h>]],
    [[return omp_get_num_threads();]]
  )])
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD \
    && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD \
    && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

# If ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no".
if test "${ac_pkg_openmp}" = no; then
  OPENMP_CFLAGS=''
  AC_MSG_WARN([OpenMP does NOT appear to be supported on this system.])
else
  OPENMP_CFLAGS="$OPENMP_CFLAGS -DHAVE_OPENMP"
fi

AC_SUBST(OPENMP_CFLAGS)

##
## MPI
##
## Based on the autoconf.ac from the Rmpi package
## <http://www.stats.uwo.ca/faculty/yu/Rmpi/>.
##
MPI=yes

# mpi discovery by argument or environment

if test -n "${MPI_INCLUDE_PATH}" ; then
  MPI_INCLUDE_PATH="${MPI_INCLUDE_PATH}"
  echo "Setting MPI include path from MPI_INCLUDE_PATH"
fi

if test -n "${MPI_LIB_PATH}" ; then
  MPI_LIB_PATH="${MPI_LIB_PATH}"
  echo "Setting MPI lib path from MPI_LIB_PATH"
fi

if test -n "${MPI_TYPE}" ; then
  MPI_TYPE="${MPI_TYPE}"
  echo "Setting MPI type from MPI_TYPE"
fi

dnl --with-mpi-include-path
AC_ARG_WITH([mpi-include-path],
  AS_HELP_STRING([--with-mpi-include-path=DIR], [location of MPI header files]),
  [
    if test "${withval}" != "no" ; then
      MPI_INCLUDE_PATH="${withval}"
    fi
  ])

dnl --with-mpi-lib-path
AC_ARG_WITH([mpi-lib-path],
  AS_HELP_STRING([--with-mpi-lib-path=DIR], [location of MPI library files]),
  [
    if test "${withval}" != "no" ; then
      MPI_LIB_PATH="${withval}"
    fi
  ])

dnl --with-mpi-type
AC_ARG_WITH([mpi-type],
  AS_HELP_STRING([--with-mpi-type=MPI_TYPE], [the type of MPI: OPENMPI, LAM or MPICH]),
  [
    if test "${withval}" != "no" ; then
      MPI_TYPE="${withval}"
    fi
  ])

if test -n "${MPI_INCLUDE_PATH}" -a -n "${MPI_LIB_PATH}" ; then
  if test ! -f "${MPI_INCLUDE_PATH}/mpi.h"; then
    AC_MSG_ERROR([Value of MPI_INCLUDE_PATH does not contain mpi.h])
  fi
  MPI_INCLUDES="-I${MPI_INCLUDE_PATH}"
  MPI_LIBS="-L${MPI_LIB_PATH}"
elif test -n "${MPI_INCLUDE_PATH}" -o -n "${MPI_LIB_PATH}" ; then
  AC_MSG_ERROR([Must specify both the include and lib path])
fi

# auto configure if not set by user
if test -z "${MPI_INCLUDE_PATH}" ; then

  dnl --with-mpi-root
  AC_ARG_WITH([mpi-root],
    AS_HELP_STRING([--with-mpi-root=DIR], [location of top-level MPI directory]),
    [
      if test "${withval}" != "no" ; then
        MPI_ROOT="${withval}"
        MPI_INCLUDE_PATH="${MPI_ROOT}/include"
        MPI_LIB_PATH="${MPI_ROOT}/lib"
      fi
    ])

  if test -z "$MPI_ROOT" ; then
    # auto configure MPI_ROOT, (MPI_TYPE, MPI_INCLUDE_PATH)
    AC_MSG_CHECKING([for mpi root directory])
    for d in /opt/lib /usr/lib /usr /usr/local/lib /usr/local \
      /usr/lib64/mpi/gcc/openmpi; do
      if test -f $d/include/mpi.h && test -d $d/lib/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        break
      elif test -f $d/include/openmpi/mpi.h && test -d $d/lib/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/include/openmpi/mpi.h && test -d $d/lib64/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/openmpi/include/mpi.h && test -d $d/openmpi/lib; then
        MPI_ROOT=$d/openmpi
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/include/mpi.h && test -d $d/lib64; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include
        break
      elif test -f $d/lam/include/mpi.h && test -f $d/lam/lib/libmpi.so; then
        MPI_ROOT=$d/lam
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/include/lam/mpi.h && test -f $d/lib/lam/libmpi.so; then
        MPI_ROOT=$d
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/include/lam/mpi.h && test -f $d/lib64/lam/libmpi.so; then
        MPI_ROOT=$d
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/mpich/include/mpi.h; then
        MPI_ROOT=$d/mpich
        MPI_TYPE="MPICH"
        MPI_INCLUDE_PATH=$d/mpich/include
        break
      elif test -f $d/mpi/include/mpi.h; then
        MPI_ROOT=$d/mpi
        MPI_INCLUDE_PATH=$d/mpi/include
        break
      elif test -f $d/include/mpi.h; then
        MPI_ROOT=$d
        break
      fi
    done
    AC_MSG_RESULT([$MPI_ROOT])
  fi

  if test -n "$MPI_ROOT"; then
    # auto configure MPI_INCLUDE_PATH
    AC_MSG_CHECKING([for mpi.h])
    for d in include include/openmpi include/lam ; do
      if test -f ${MPI_ROOT}/$d/mpi.h ; then
        AC_MSG_RESULT([${MPI_ROOT}/$d])
        MPI_INCLUDE_PATH="${MPI_ROOT}/$d"
        break
      fi
    done
    if test -z "$MPI_INCLUDE_PATH" ; then
      AC_MSG_RESULT([no])
      MPI=no
    fi
  else
    AC_CHECK_HEADER(mpi.h, , [MPI=no])
  fi
  if test -n "$MPI_INCLUDE_PATH" ; then
    MPI_INCLUDES="-I${MPI_INCLUDE_PATH}"
  fi

  dnl auto configure MPI_LIB_PATH, (MPI_INCLUDES)
  AC_MSG_CHECKING([for libmpi.{a,so,dylib}])
  if test -f ${MPI_ROOT}/lib/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  elif test -f ${MPI_ROOT}/libmpi.so; then
    dnl Fedora 10
    MPI_LIB_PATH="${MPI_ROOT}"
  elif test -f ${MPI_ROOT}/lib/openmpi/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib/openmpi"
  elif test -f ${MPI_ROOT}/lib64/openmpi/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64/openmpi"
  elif test -f ${MPI_ROOT}/lib64/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64"
  elif test -f ${MPI_ROOT}/lib/lam/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib/lam"
    MPI_INCLUDES="$MPI_INCLUDES -I${MPI_INCLUDE_PATH}/32"
  elif test -f ${MPI_ROOT}/lib64/lam/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64/lam"
    MPI_INCLUDES="$MPI_INCLUDES -I${MPI_INCLUDE_PATH}/64"
  elif test -f ${MPI_ROOT}/lib/libmpich.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib"
    AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],,[#include <mpi.h>])
  elif test -f ${MPI_ROOT}/lib64/libmpich.a; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64"
    AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],,[#include <mpi.h>])
  elif test -f ${MPI_ROOT}/lib/libmpi.dylib; then
    dnl OSX (10.6)
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  else
    AC_CHECK_LIB(mpi, main, [MPI_LIB_PATH="${MPI_ROOT}/lib"], )
  fi
  if test -n "$MPI_LIB_PATH" ; then
    AC_MSG_RESULT([$MPI_LIB_PATH])
  fi
fi

# end auto configure type
if test -z $MPI_TYPE ; then
  if   test -d ${MPI_ROOT}/openmpi;          then MPI_TYPE="OPENMPI"
  elif test -d ${MPI_ROOT}/lib/openmpi;      then MPI_TYPE="OPENMPI"
  elif test -d ${MPI_ROOT}/lib64/openmpi;    then MPI_TYPE="OPENMPI"
  elif test -f ${MPI_ROOT}/lib/liblam.so;    then MPI_TYPE="LAM"
  elif test -d ${MPI_ROOT}/lib/lam;          then MPI_TYPE="LAM"
  elif test -d ${MPI_ROOT}/lib64/lam;        then MPI_TYPE="LAM"
  elif test -f ${MPI_ROOT}/lib/libmpich.a;   then MPI_TYPE="MPICH"
  elif test -f ${MPI_ROOT}/lib64/libmpich.a; then MPI_TYPE="MPICH"
  fi
fi

if test "$MPI_TYPE" = "OPENMPI" ; then
  AC_CHECK_PROG(ORTED, orted, yes, no)
  if test "$ORTED" = no ; then
    AC_MSG_NOTICE([Cannot find orted. RSiena needs orted to run.])
    MPI=no
  fi
fi

if test "$MPI_TYPE" = "LAM" ; then
  echo "Try to find liblam.so ..."
  if test -f ${MPI_ROOT}/lib/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib"
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  elif test -f ${MPI_ROOT}/lib/lam/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib/lam"
    MPI_LIB_PATH="${MPI_ROOT}/lib/lam"
  elif test -f ${MPI_ROOT}/lib64/lam/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib64/lam"
    MPI_LIB_PATH="${MPI_ROOT}/lib64/lam"
  else
    AC_CHECK_LIB(lam, main, [ echo "liblam not found. Probably not LAM-MPI" ])
  fi
fi

if test -n "$MPI_LIB_PATH" ; then
  MPI_LIBS="-L${MPI_LIB_PATH}"
else
  MPI=no
fi

AC_CHECK_LIB(util, openpty, [ MPI_LIBS="$MPI_LIBS -lutil" ])
AC_CHECK_LIB(pthread, main, [ MPI_LIBS="$MPI_LIBS -lpthread" ])

MPI_DEFS="-D$MPI_TYPE"
if test "$MPI_TYPE" = "OPENMPI" -o "$MPI_TYPE" = "LAM" -o "$MPI_TYPE" = "MPICH2"; then
  MPI_DEFS="$MPI_DEFS -DMPI2"
fi

OSTYPE="`uname`"
if test "$OSTYPE" = "Darwin"; then
  MPI_DEFS="$MPI_DEFS -DMAC"
fi

if test "$MPI" = "yes" ; then

PKG_CPPFLAGS="$MPI_INCLUDES $MPI_DEFS"

  case "$MPI_TYPE" in
    OPENMPI) PKG_LIBS="-lmpi_cxx -lmpi ${MPI_LIBS}";;
    LAM)     PKG_LIBS="-lmpi_cxx -llam ${MPI_LIBS}";;
    MPICH)   PKG_LIBS="-lmpich -lmpl -lopa ${MPI_LIBS}";;
    MPICH2)  PKG_LIBS="-lmpich -lmpl -lopa ${MPI_LIBS}";;
    *)
      AC_MSG_NOTICE([Unknown type of mpi, compiling without mpi support. Use --with-mpi-xxxx to specify it.])
      PKG_CPPFLAGS=
      PKG_LIBS=
      ;;
  esac
else
  PKG_CPPFLAGS=
  PKG_LIBS=
fi

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

# substitute externalized source list
PKG_SOURCES=$(cat src/sources.list)
AC_SUBST(PKG_SOURCES)

AC_OUTPUT(src/Makevars)
