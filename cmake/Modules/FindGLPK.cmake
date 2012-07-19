find_package(PkgConfig)
pkg_check_modules(PC_GLPK glpk)
set(GLPK_DEFINITIONS ${PC_GLPK_CFLAGS_OTHER})

FIND_PATH(GLPK_INCLUDE_DIR glpk.h
    HINTS 
        ${PC_GLPK_INCLUDEDIR} 
        ${PC_GLPK_INCLUDE_DIRS} 
        $ENV{GLPK_DIR}/include 
        $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES glpk)
set(GLPK_INCLUDE_DIRS ${GLPK_INCLUDE_DIR})

FIND_LIBRARY(GLPK_LIBRARY glpk 
    HINTS
        ${PC_GLPK_LIBDIR} 
        ${PC_GLPK_LIBRARY_DIRS} 
        $ENV{GLPK_DIR}/lib)
set(GLPK_LIBRARIES ${GLPK_LIBRARY})

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set GLPK_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(GLPK DEFAULT_MSG
                                  GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)
