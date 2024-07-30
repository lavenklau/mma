cmake_minimum_required(VERSION 3.8)

project(mma)

enable_language(CUDA)
# disable all warning for clean output
add_definitions(-w)
find_package(CUDA REQUIRED)

# Do what the new package does
find_library(CUDA_DRIVER_LIBRARY
             NAMES cuda_driver cuda
             HINTS ${CUDA_TOOLKIT_ROOT_DIR}
                   ENV CUDA_PATH
             PATH_SUFFIXES nvidia/current lib64 lib/x64 lib)
if (NOT CUDA_DRIVER_LIBRARY)
    # Don't try any stub directories until we have exhausted all other search locations.
    find_library(CUDA_DRIVER_LIBRARY
                 NAMES cuda_driver cuda
                 HINTS ${CUDA_TOOLKIT_ROOT_DIR}
                       ENV CUDA_PATH
                 PATH_SUFFIXES lib64/stubs lib/x64/stubs lib/stubs stubs)
endif ()

find_package(Eigen3 REQUIRED)


set (BUILD_TEST ON CACHE BOOL "whether to build test target")
set (USE_MATLAB OFF CACHE BOOL "whether to use matlab support, usefull when debuging")

if(BUILD_TEST)
    # set (USE_MATLAB ON CACHE BOOL STRING FORCE)
endif()

set (CMAKE_CXX_STANDARD 17)

set (CMAKE_VERBOSE_MAKEFILE ON)

list (APPEND libfiles mmaOpt.cpp mmaOpt.cu lib.cu)

if (USE_MATLAB)
    find_package(Matlab REQUIRED COMPONENTS ENG_LIBRARY MX_LIBRARY MAT_LIBRARY)
    add_definitions(-DENABLE_MATLAB)
    add_definitions(-DGVECTOR_WITH_MATLAB)
    message(STATUS "matlab eng library ${Matlab_ENG_LIBRARY}")
    message(STATUS "matlab mx library ${Matlab_MX_LIBRARY}")
    list(APPEND libfiles matlab_utils.cpp gpuVector.cpp)
endif()

add_library(mma SHARED ${libfiles})
set_target_properties(mma PROPERTIES CUDA_ARCHITECTURES "61;75;86")
target_include_directories(mma PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_compile_options(mma PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    --extended-lambda
    --std=c++17
>)

if(USE_MATLAB)
    target_include_directories(mma PRIVATE ${Matlab_INCLUDE_DIRS})
    target_link_libraries(mma PRIVATE ${Matlab_MX_LIBRARY})
    target_link_libraries(mma PRIVATE ${Matlab_ENG_LIBRARY})
    target_link_libraries(mma PRIVATE ${Matlab_MAT_LIBRARY})
endif()
target_include_directories(mma PRIVATE "${CMAKE_SOURCE_DIR}/Common")
target_compile_definitions(mma PRIVATE -D_EXPORT_MMAOPT)
target_link_libraries (mma PRIVATE Eigen3::Eigen)
target_link_libraries(mma PRIVATE ${CUDA_LIBRARIES})
target_link_libraries(mma PRIVATE ${CUDA_DRIVER_LIBRARY})
target_link_libraries(mma PRIVATE ${CUDA_cusolver_LIBRARY})

if(BUILD_TEST)
    if(USE_MATLAB)
        add_executable(test_mma test_suits.cpp matlab_utils.cpp)
    else()
        add_executable(test_mma test_suits.cpp)
    endif()
    target_include_directories(test_mma PRIVATE ${Matlab_INCLUDE_DIRS})
    target_link_libraries(test_mma PRIVATE mma)
    target_link_libraries(test_mma PRIVATE Eigen3::Eigen)
    target_link_libraries(test_mma PRIVATE ${Matlab_MX_LIBRARY})
    target_link_libraries(test_mma PRIVATE ${Matlab_ENG_LIBRARY})
    target_link_libraries(test_mma PRIVATE ${Matlab_MAT_LIBRARY})
endif()

