dnl ### configure.ac					-*- Autoconf -*-

dnl dnl Based on OpenMP detection in R's configure script,
dnl dnl which in turn is based on autoconf's openmp.m4, 
dnl dnl with is licensed under GPL-3 with the 
dnl dnl Autoconf Configure Script Exception, version 3.0.

AC_PREREQ(2.62)

dnl dnl derived from the BART R package v 2.3
AC_INIT([BART],[2.3],[],[],[])

dnl CXX11FLAGS=`R CMD config CXX11FLAGS`
dnl CXX11=`R CMD config CXX11`
dnl CXX11STD=`R CMD config CXX11STD`
dnl LDFLAGS=`R CMD config LDFLAGS`

dnl CXX="${CXX11}"
dnl CXXFLAGS="${CXXFLAGS} ${CXX11STD} ${CXX11FLAGS}"
dnl CPPFLAGS="${CPPFLAGS} -D_OPENMP"

dnl AC_LANG(C++)
dnl dnl the meat of R's m4/openmp.m4
dnl   OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
dnl   AC_ARG_ENABLE([openmp],
dnl     [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
dnl   if test "$enable_openmp" != no; then
dnl     AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
dnl       [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
dnl       [AC_LINK_IFELSE([_AC_LANG_OPENMP],
dnl 	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
dnl 	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
dnl 	  for ac_option in -xopenmp -fopenmp -qopenmp \
dnl                            -openmp -mp -omp -qsmp=omp -homp \
dnl 			   -fopenmp=libomp \
dnl                            -Popenmp --openmp \
dnl 			   "-Xclang -fopenmp -Xlinker -lomp"; do
dnl 	    ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
dnl 	    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
dnl 	    AC_LINK_IFELSE([_AC_LANG_OPENMP],
dnl 	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
dnl 	    _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
dnl 	    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
dnl 	      break
dnl 	    fi
dnl 	  done])])
dnl     case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
dnl       "none needed" | unsupported)
dnl 	;; #(
dnl       *)
dnl 	OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
dnl     esac
dnl   fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

