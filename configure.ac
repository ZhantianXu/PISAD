AC_INIT([PISAD], [1.0], [234712280@csu.edu.cn])
AM_INIT_AUTOMAKE([-Wall -Werror foreign no-dependencies])

AC_CONFIG_HEADERS([config.h])

# check gcc
AC_PROG_CXX
AC_PROG_CC
AM_PROG_AR

# check Boost 
#AC_CHECK_LIB([boost_system], [boost::system])

# check HDF5 
#AC_CHECK_LIB([hdf5], [H5Fopen])

# check tcmalloc
#AC_CHECK_LIB([tcmalloc], [malloc])

# check openmp
#AX_OPENMP

CXXFLAGS="-g -gdwarf-4 -O3"

# generate Makefile
AC_CONFIG_FILES([Makefile src/Makefile snpcalling/Makefile])
AC_OUTPUT
