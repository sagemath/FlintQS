AC_PREREQ([2.66])
AC_INIT([FlintQS], [1.0], [sage-devel@googlegroups.com])

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])

AC_CONFIG_FILES([
    Makefile
    src/Makefile
])


AC_CHECK_HEADERS([gmp.h])
AC_CHECK_LIB([gmp], [__gmpz_init], ,
    AC_MSG_ERROR([GNU Multiple Precision Arithmetic Library not found.
        Set CFLAGS/LDFLAGS if it is installed in a non-standard directory])
)

AC_PROG_CXX
AC_OUTPUT
