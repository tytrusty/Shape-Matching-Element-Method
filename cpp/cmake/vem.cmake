#this file setups up targets
#can be imported into other projects to use vem

# vem Options
set(vem_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(vem_SOURCE_DIR "${vem_ROOT}")
set(vem_INCLUDE_DIR ${vem_ROOT}/include)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${vem_SOURCE_DIR}/cmake)

option(VEM_USE_OPENMP OFF)
if(NOT APPLE)
    option(LIBIGL_WITH_EMBREE ON)
else ()
    option(LIBIGL_WITH_EMBREE "Use Embree" ON)
endif()

find_package(LIBIGL REQUIRED)

if(VEM_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

add_library(vem INTERFACE)

include_directories(${vem_INCLUDE_DIR})

target_include_directories(vem INTERFACE ${vem_INCLUDE_DIR})
target_compile_definitions(vem INTERFACE -DSIM_DATA_DIRECTORY=${vem_ROOT}/models)
target_link_libraries(vem INTERFACE igl::core;igl::embree)

if(OpenMP_CXX_FOUND)

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        execute_process(
            COMMAND brew --prefix libomp 
            RESULT_VARIABLE BREW_OMP
            OUTPUT_VARIABLE BREW_OMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        include_directories(${BREW_OMP_PREFIX}/include)

        target_link_libraries(vem INTERFACE OpenMP::OpenMP_CXX)
        target_compile_definitions(vem INTERFACE -DVEM_USE_OPENMP)

    elseif()
        include_directories(${OpenMP_CXX_INCLUDE_DIRS})

        target_link_libraries(vem INTERFACE OpenMP::OpenMP_CXX)
        target_compile_definitions(vem INTERFACE -DVEM_USE_OPENMP)
    endif()
    
endif()