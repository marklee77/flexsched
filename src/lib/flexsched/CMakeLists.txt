#find_package(GLPK REQUIRED)

#include_directories(${GLPK_INCLUDE_DIR})

add_library(flexsched 
            util.c 
            greedy/greedy.c 
            vector/vecutil.c vector/hetero.c)

#target_link_libraries(flexsched ${GLPK_LIBRARIES})
