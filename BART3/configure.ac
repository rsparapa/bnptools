### configure.ac					-*- Autoconf -*-

dnl Based on OpenMP detection in R's configure script,
dnl which in turn is based on autoconf's openmp.m4, 
dnl with is licensed under GPL-3 with the 
dnl Autoconf Configure Script Exception, version 3.0.

AC_PREREQ(2.62)

dnl derived from the BART R package v 2.3
AC_INIT([BART],[2.3],[],[],[])

dnl CXX11FLAGS=`R CMD config CXX11FLAGS`
dnl CXX11=`R CMD config CXX11`
dnl CXX11STD=`R CMD config CXX11STD`
LDFLAGS=`R CMD config LDFLAGS`

dnl CXX="${CXX11}"
dnl CXXFLAGS="${CXXFLAGS} ${CXX11STD} ${CXX11FLAGS}"
CPPFLAGS="${CPPFLAGS} -D_OPENMP"

AC_LANG(C++)
dnl the meat of R's m4/openmp.m4
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
	  for ac_option in -xopenmp -fopenmp -qopenmp \
                           -openmp -mp -omp -qsmp=omp -homp \
			   -fopenmp=libomp \
                           -Popenmp --openmp \
			   "-Xclang -fopenmp -Xlinker -lomp"; do
	    ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
	    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
	    AC_LINK_IFELSE([_AC_LANG_OPENMP],
	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
	    _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
	    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
	      break
	    fi
	  done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
	;; #(
      *)
	OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
  fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

