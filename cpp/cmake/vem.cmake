#this file setups up targets
#can be imported into other projects to use vem

# vem Options
set(vem_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(vem_SOURCE_DIR "${vem_ROOT}")
set(vem_INCLUDE_DIR ${vem_ROOT}/include)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${vem_SOURCE_DIR}/cmake)

option(VEM_USE_OPENMP OFF)

find_package(LIBIGL REQUIRED)

if(VEM_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# Dependencies are linked as INTERFACE targets unless libigl is compiled as a static library
if(vem_USE_STATIC_LIBRARY)
    
    file(GLOB SOURCES_vem "${vem_SOURCE_DIR}/src/*.cpp")

    add_library(vem STATIC ${SOURCES_vem})
    
    if(MSVC)
        target_compile_options(vem PRIVATE /w) # disable all warnings (not ideal but...)
    else()
        #target_compile_options(${module_libname} PRIVATE -w) # disable all warnings (not ideal but...)
    endif()
else()
    add_library(vem INTERFACE)
endif()

include_directories(${vem_INCLUDE_DIR})

#prepreprocessor definition for static library 
if(vem_USE_STATIC_LIBRARY)
    target_include_directories(vem PUBLIC ${vem_INCLUDE_DIR})
    target_compile_definitions(vem PUBLIC -DSIM_DATA_DIRECTORY=${vem_ROOT}/models)
    target_link_libraries(vem PUBLIC igl::core)

    target_compile_definitions(vem PUBLIC -DSIM_STATIC_LIBRARY)
else()
    target_include_directories(vem INTERFACE ${vem_INCLUDE_DIR})
    target_compile_definitions(vem INTERFACE -DSIM_DATA_DIRECTORY=${vem_ROOT}/models)
    target_link_libraries(vem INTERFACE igl::core)
endif()

if(OpenMP_CXX_FOUND)

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        execute_process(
            COMMAND brew --prefix libomp 
            RESULT_VARIABLE BREW_OMP
            OUTPUT_VARIABLE BREW_OMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        include_directories(${BREW_OMP_PREFIX}/include)

        if(vem_USE_STATIC_LIBRARY)
            target_link_libraries(vem PUBLIC OpenMP::OpenMP_CXX)
            target_compile_definitions(vem PUBLIC -DVEM_USE_OPENMP)
        else()
            target_link_libraries(vem INTERFACE OpenMP::OpenMP_CXX)
            target_compile_definitions(vem INTERFACE -DVEM_USE_OPENMP)
        endif()

    elseif()
        include_directories(${OpenMP_CXX_INCLUDE_DIRS})

        if(vem_USE_STATIC_LIBRARY)
            target_link_libraries(vem PUBLIC OpenMP::OpenMP_CXX)
            target_compile_definitions(vem PUBLIC -DVEM_USE_OPENMP)
        else()
            target_link_libraries(vem INTERFACE OpenMP::OpenMP_CXX)
            target_compile_definitions(vem INTERFACE -DVEM_USE_OPENMP)
        endif()
    endif()
    
endif()