FIND_PATH(GLPK_INCLUDE_DIR
    NAMES glpk.h
    PATHS /usr/include /usr/include/glpk $ENV{GLPK_DIR}/include
    DOC "Directory where GLPK header files are stored")

FIND_LIBRARY(GLPK_LIBRARY
    NAMES glpk 
    PATHS /usr/lib $ENV{GLPK_DIR}/lib)

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set GLPK_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(GLPK DEFAULT_MSG
                                  GLPK_LIBRARY GLPK_INCLUDE_DIR)

mark_as_advanced(GLPK_INCLUDE_DIR GLPK_LIBRARY)
