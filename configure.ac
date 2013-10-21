AC_INIT([hglmCPP], 0.01)
AC_CONFIG_AUX_DIR([tools])

AC_ARG_WITH([compiler-bindir],
	AC_HELP_STRING([--with-compiler-bindir=compiler_bindir],
		[specify the directory containing gcc, e.g. /usr/bin]),
	[COMPILER_BINDIR="$withval"])

AC_ARG_WITH([lapack],
	AC_HELP_STRING([--with-lapack=lapack],
		[specify your LAPACK library, e.g. openblas, atlas or mkl]),
	[LAPACK="$withval"],[LAPACK="internal"])
	
AC_ARG_WITH([lapack-include],
	AC_HELP_STRING([--with-lapack-include=lapack_include],
		[specify the LAPACK include directory, e.g. -I/usr/local/openblas/include]),
	[LAPACK_INCLUDE="$withval"])
	
AC_ARG_WITH([lapack-flags],
	AC_HELP_STRING([--with-lapack-flags=lapack_flags],
		[specify the LAPACK compiler flags, e.g. -m64]),
	[LAPACK_FLAGS="$withval"])
		
AC_ARG_WITH([lapack-library],
	AC_HELP_STRING([--with-lapack-library=lapack_library],
		[specify the LAPACK library directory, e.g. -L/usr/local/openblas/lib]),
	[LAPACK_LIBRARY="$withval"])

AC_ARG_WITH([lapack-options],
	AC_HELP_STRING([--with-lapack-options=lapack_options],
		[specify the LAPACK link options, e.g. -lopenblas]),
	[LAPACK_OPTIONS="$withval"])

AC_ARG_WITH([r-include],
	AC_HELP_STRING([--with-r-include=r_include],
		[specify the full path to your R header files, e.g. /usr/share/R/include]),
	[R_INCLUDE="$withval"])

AC_ARG_WITH([r-lib],
	AC_HELP_STRING([--with-r-lib=r-lib],
		[specify the full path to your R shared libraries, e.g. /usr/lib/R/lib]),
	[R_LIB="$withval"])

AC_MSG_RESULT([---------------------------])
AC_MSG_RESULT([Checking for critical tools])
AC_MSG_RESULT([---------------------------])

################################
## Testing existence of tools ##
################################

## Check host
AC_CANONICAL_HOST

## Check if gcc exists
AC_MSG_CHECKING(for environment variable COMPILER_BINDIR)
if test -z "${COMPILER_BINDIR}"; then
	AC_MSG_RESULT(ok)
else
	AC_MSG_RESULT(ok)
	AC_CHECK_FILE([${COMPILER_BINDIR}/gcc], [HAS_GCC="yes"])
	if test -z "${HAS_GCC}"; then
		AC_MSG_ERROR(failed)
		MISC_PARAMS=""
	else
		AC_MSG_RESULT(ok)
		MISC_PARAMS="--compiler-bindir ${COMPILER_BINDIR}"
	fi
fi

## Check where R is installed
AC_MSG_CHECKING(for R install path)
if test -z "${R_HOME}"; then
	AC_MSG_ERROR(failed)
else
	AC_MSG_RESULT(ok)
	CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
	AC_MSG_CHECKING(for environment variable R_INCLUDE)
	if test -z "${R_INCLUDE}"; then
		R_INCLUDE=`"${R_HOME}/bin/R" CMD config --cppflags`
		AC_MSG_RESULT(ok)
	else
		R_INCLUDE="-I${R_INCLUDE}"
		AC_MSG_RESULT(ok)
	fi
	AC_MSG_CHECKING(for environment variable R_LIB)
	if test -z "${R_LIB}"; then
		R_LIB=`"${R_HOME}/bin/R" CMD config --ldflags`
		AC_MSG_RESULT(ok)
	else
		R_LIB="-L${R_LIB} -lR"
		AC_MSG_RESULT(ok)
	fi
fi

#########################################
## Make variable available in makefile ##
#########################################

AC_SUBST(LAPACK)
AC_SUBST(LAPACK_INCLUDE)
AC_SUBST(LAPACK_FLAGS)
AC_SUBST(LAPACK_LIBRARY)
AC_SUBST(LAPACK_OPTIONS)
AC_SUBST(MISC_PARAMS)
AC_SUBST(R_HOME)
AC_SUBST(R_INCLUDE)
AC_SUBST(R_LIB)
AC_SUBST(CPICFLAGS)

AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT

AC_MSG_RESULT([-------------------------])
AC_MSG_RESULT([Installing FMPR R-package])
AC_MSG_RESULT([-------------------------])
