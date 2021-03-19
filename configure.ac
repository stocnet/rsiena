#
# After changing this file run `autoconf` to produce the 'configure' script. 
#
# - Run with automatic detection: `R CMD INSTALL {pkg}`
#

AC_PREREQ([2.69])
AC_INIT([RSiena], 1.2.32)
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

AC_OUTPUT(src/Makevars)

echo "
 --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}
  ================
  
    cppflags: ${CPPFLAGS} 
    cxxflags: ${CXXFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
