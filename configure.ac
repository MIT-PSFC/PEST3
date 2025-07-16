###################################################################
# $Id: $
###################################################################

AC_INIT(pest3code,3.3,pletzer@txcorp.com)
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

dnl ######################################################################
dnl
dnl Compilers
dnl
dnl ######################################################################

AC_PROG_FC
AC_PROG_CXX
AC_PROG_CC

LT_INIT
AC_PROG_LIBTOOL

AC_PROG_INSTALL

# Default prefix
AC_PREFIX_DEFAULT(/internal/nubeam)

# Do not automatically regenerate Makefile.in's
AM_MAINTAINER_MODE

# Prevent testing executables as fails on parallel machines

# all subsequent test files will be .f90
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT([f90])

dnl###################################################################
dnl BLAS
dnl###################################################################

searchdirs="/usr/lib /usr/lib64"
BLAS_LIB="NotSet"
AC_ARG_WITH([blas-lib],
        [AS_HELP_STRING([--with-blas-lib=<blas-lib>],
        [full path to BLAS library])],
	BLAS_LIB="$withval")
if test "$BLAS_LIB" = "NotSet"; then
  for dir in $searchdirs; do	
     AC_CHECK_FILE($dir/libblas.a, 
     BLAS_LIB="$dir/libblas.a")	
  done
fi
if test "$BLAS_LIB" = "NotSet"; then
  AC_MSG_ERROR([Cannot find location of BLAS library. Consider using --with-blas-lib=<...>])
fi

AC_SUBST([BLAS_LIB],["$BLAS_LIB"])

dnl###################################################################
dnl LAPACK
dnl###################################################################

searchdirs="/usr/lib /usr/lib64"
LAPACK_LIB="NotSet"
AC_ARG_WITH([lapack-lib],
        [AS_HELP_STRING([--with-lapack-lib=<lapack-lib>],
        [full path to LAPACK library])],
	LAPACK_LIB="$withval")
if test "$LAPACK_LIB" = "NotSet"; then
  for dir in $searchdirs; do	
     AC_CHECK_FILE($dir/liblapack.a, 
     LAPACK_LIB="$dir/liblapack.a")	
  done
fi
if test "$LAPACK_LIB" = "NotSet"; then
  AC_MSG_ERROR([Cannot find location of LAPACK library. Consider using --with-lapack-lib=<...>])
fi

AC_SUBST([LAPACK_LIB],["$LAPACK_LIB"])


dnl###################################################################
dnl NetCDF
dnl###################################################################

searchdirs="/usr/local/ /usr $HOME/software /opt/netcdf"
NETCDF_DIR="NotSet"
AC_ARG_WITH([netcdf-dir],
        [AS_HELP_STRING([--with-netcdf-dir=<netcdf-dir>],
        [location of NETCDF directory])],
	NETCDF_DIR="$withval")
if test "$NETCDF_DIR" = "NotSet"; then
  for dir in $searchdirs; do	
     AC_CHECK_FILE($dir/lib/libpspline.a, 
     NETCDF_DIR="$dir")	
  done
fi
if test "$NETCDF_DIR" = "NotSet"; then
  AC_MSG_ERROR([Cannot find location of NETCDF libraries. Consider using --with-netcdf-dir=<...>])
fi

AC_SUBST([NETCDF_DIR],["$NETCDF_DIR"])

dnl ##################################################################
dnl MDSPLUS libraries
dnl ##################################################################

searchdirs="/usr/mdsplus /usr/local/mdsplus /opt/mdsplus"
MDSPLUS_LIBS=""
MDSPLUS_DIR="NotSet"
AC_ARG_WITH([mdsplus-dir],
        [AS_HELP_STRING([--with-mdsplus-dir=<mdsplus-dir>],
        [location of MDSPLUS directory])],
	MDSPLUS_DIR="$withval")
if test "$MDSPLUS_DIR" = "NotSet"; then
  for dir in $searchdirs; do	
     AC_CHECK_FILES($dir/lib/libTreeShr.so, 
     MDSPLUS_DIR="$dir")	
  done
fi
if test "$MDSPLUS_DIR" = "NotSet"; then
  AC_MSG_WARN([Could not find location of MDSPLUS libraries. Consider using --with-mdsplus-dir=<...>])
else
  MDSPLUS_LIBS="-L$MDSPLUS_DIR -lTreeShr -lMdsShr -lTdiShr"
fi

AC_SUBST([MDSPLUS_LIBS],["$MDSPLUS_LIBS"])

dnl ##################################################################
dnl NTCC libraries
dnl ##################################################################

searchdirs="$HOME/ntcc /usr/local/ntcc"
NTCC_DIR="NotSet"
AC_ARG_WITH([ntcc-dir],
        [AS_HELP_STRING([--with-ntcc-dir=<ntcc-dir>],
        [location of NTCC directory])],
	NTCC_DIR="$withval")
if test "$NTCC_DIR" = "NotSet"; then
  for dir in $searchdirs; do	
     AC_CHECK_FILE($dir/lib/libpspline.a, 
     NTCC_DIR="$dir")	
  done
fi
if test "$NTCC_DIR" = "NotSet"; then
  AC_MSG_ERROR([Cannot find location of NTCC libraries. Consider using --with-ntcc-dir=<...>])
fi

AC_SUBST([NTCC_DIR],["$NTCC_DIR"])

######################################################################
# Mixed Fortran C/C++
######################################################################
# determine the fortran libraries
AC_FC_LIBRARY_LDFLAGS

# determine the name of the object file that initializes fortran
# IO, etc
#builtin(include, config/tx_fc_find_main.m4)
# this will set $TX_FC_MAIN_OBJS
#TX_FC_FIND_MAIN_OBJS

######################################################################
# Create Makefiles
######################################################################
AC_CONFIG_FILES(
[Makefile]
[pest3/Makefile]
)


AC_OUTPUT
