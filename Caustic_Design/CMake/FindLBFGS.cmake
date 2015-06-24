# Try to find LBFGS
# Once done this will define
#  LBFGS_FOUND        - Indication as to whether LBFGS was found
#  LBFGS_INCLUDE_DIRS - The LBFGS include directories
#  LBFGS_LIBRARIES    - The libraries needed to use LBFGS
#  LBFGS_DEFNINITIONS - Compiler flags needed to use LBFGS

find_package (PkgConfig)
pkg_check_modules (PC_LBFGS QUIET liblbfgs)
set (LBFGS_DEFINITIONS ${PC_LBFGS_CFLAGS_OTHER})

find_path(LBFGS_INCLUDE_DIR NAMES lbfgs.h
    HINTS @CMAKE_INSTALL_PREFIX@/include ${PC_LBFGS_INCLUDEDIR} ${PC_LBFGS_INCLUDE_DIRS} /usr/include /usr/local/include /opt/local/include)

find_library (LBFGS_LIBRARY NAMES liblbfgs lbfgs
    HINTS @CMAKE_INSTALL_PREFIX@/lib ${PC_LBFGS_LIBDIR} ${PC_LBFGS_LIBRARY_DIRS} /usr/lib /usr/local/lib /opt/local/lib)

set (LBFGS_LIBRARIES ${LBFGS_LIBRARY})
set (LBFGS_INCLUDE_DIRS ${LBFGS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args (LBFGS DEFAULT_MSG
    LBFGS_LIBRARY LBFGS_INCLUDE_DIR)

mark_as_advanced (LBFGS_INCLUDE_DIR LBFGS_LIBRARY)