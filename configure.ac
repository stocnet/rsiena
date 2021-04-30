#
# After changing this file run `autoconf` to produce the 'configure' script. 
#
# - Run with automatic detection: `R CMD INSTALL {pkg}`
#

AC_PREREQ([2.69])
AC_INIT([RSiena], @VERSION@)
AC_CONFIG_SRCDIR([src])

# Find the compiler and compiler flags used by R.
R_HOME=`R RHOME`
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi
RBIN="${R_HOME}/bin/R"
LIBS="${PKG_LIBS}"

# Checking for C++
# base on: https://unconj.ca/blog/an-autoconf-primer-for-r-package-authors.html 
CXX=`"${RBIN}" CMD config CXX`
CXXFLAGS=`"${RBIN}" CMD config CXXFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# Check for typedefs, structures, and compiler characteristics
AC_CHECK_HEADERS([sys/time.h])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

# substitute externalized source list
PKG_SOURCES=$(cat src/sources.list)
AC_SUBST(PKG_SOURCES)

# copy from GitHub repository neurodata/knorR repository
AC_CONFIG_FILES([src/Makevars.tmp:src/Makevars.in], [
    if test -f src/Makevars && cmp -s src/Makevars.tmp src/Makevars; then
      AC_MSG_NOTICE([creating src/Makevars])
      AC_MSG_NOTICE([src/Makevars is unchanged])
      rm src/Makevars.tmp
    else
      AC_MSG_NOTICE([creating src/Makevars])
      mv src/Makevars.tmp src/Makevars
    fi
  ]
)

AC_OUTPUT

echo "
 --------------------------------------------------
  Configuration for ${PACKAGE_NAME}
  ================
  
    cppflags: ${CPPFLAGS} 
    cxxflags: ${CXXFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
