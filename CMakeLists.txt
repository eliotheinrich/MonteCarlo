cmake_minimum_required(VERSION 3.15..3.26)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(MonteCarlo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message("Build type: " ${CMAKE_BUILD_TYPE})
message("Compiler: " ${CMAKE_CXX_COMPILER})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-deprecated-declarations -pthread -fopenmp")

# Ignore annoying Intel warnings
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable=2196")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
string(REPLACE ":" ";" CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
list(APPEND CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})

# -------------------------------- #
# ------ Fetching externals ------ #
# -------------------------------- #

foreach(PATH IN LISTS CMAKE_PREFIX_PATH)
    include_directories("${PATH}/include")
    link_directories("${PATH}/lib")
endforeach()

add_definitions(-DFMT_HEADER_ONLY)

include(FetchContent)

set(EXTERNALS "")

# Python wrappers
find_package(Python 3.10 
    REQUIRED COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule
)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR
)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

if (DEFINED QUTILS_SOURCE) 
    set(qutils_SOURCE_DIR ${QUTILS_SOURCE})
    add_subdirectory(${qutils_SOURCE_DIR} ${CMAKE_BINARY_DIR}/qutils-build)
else()
    FetchContent_Declare(
        qutils
        GIT_REPOSITORY https://github.com/eliotheinrich/qutils.git
        GIT_TAG main
        GIT_SHALLOW TRUE
    )
    list(APPEND EXTERNALS qutils)
endif()

if (DEFINED UTILS_SOURCE) 
    set(utils_SOURCE_DIR ${UTILS_SOURCE})
    add_subdirectory(${utils_SOURCE_DIR} ${CMAKE_BINARY_DIR}/utils-build)
else()
    FetchContent_Declare(
        utils
        GIT_REPOSITORY https://github.com/eliotheinrich/utils.git
        GIT_TAG main
        GIT_SHALLOW TRUE
    )
    list(APPEND EXTERNALS utils)
endif()

if (DEFINED DATAFRAME_SOURCE)
    set(dataframe_SOURCE_DIR ${DATAFRAME_SOURCE})
else()
    FetchContent_Declare(
        dataframe
        GIT_REPOSITORY https://github.com/eliotheinrich/dataframe
        GIT_TAG master
        GIT_SHALLOW TRUE
    )
    list(APPEND EXTERNALS dataframe)
endif()

FetchContent_Declare(
  glaze
  GIT_REPOSITORY https://github.com/stephenberry/glaze.git
  GIT_TAG cc4ce0efe5133c23ec55f07d9e753f7dd1abb0b8
)
list(APPEND EXTERNALS glaze)

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG master
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
list(APPEND EXTERNALS eigen)

set(FINUFFT_FFTW_LIBRARIES "${CONDA_PREFIX}/lib/libfftw3.so" "${CONDA_PREFIX}/lib/libfftw3.a")

FetchContent_Declare(
  finufft
  GIT_REPOSITORY https://github.com/flatironinstitute/finufft.git
  GIT_TAG master
)

list(APPEND EXTERNALS finufft)

FetchContent_MakeAvailable(${EXTERNALS})

include_directories(
    ${glaze_SOURCE_DIR}/include
    ${dataframe_SOURCE_DIR}/src
    ${utils_SOURCE_DIR}/src
    ${qutils_SOURCE_DIR}/src
    ${eigen_SOURCE_DIR}
)

set(ISING_DIR src/Ising)
set(MB_ISING_DIR src/MultibodyIsing)
set(SPIN2D_DIR src/Spin2d)
set(SPIN3D_DIR src/Spin3d)
set(GRAPH_DIR src/Graph)

set(MODEL_DIRS ${ISING_DIR} ${MB_ISING_DIR} ${SPIN2D_DIR} ${SPIN3D_DIR} ${GRAPH_DIR})
foreach(DIR IN LISTS MODEL_DIRS)
	add_subdirectory(${DIR})
endforeach()

set(SIMULATOR_LIBS ising multibody_ising mc_graph spin2d spin3d)
target_link_libraries(spin3d PUBLIC finufft)

if (MC_BUILDING_PYTHON)
    nanobind_add_module(
        montecarlo_bindings
        STABLE_ABI
        NB_STATIC
        src/PyMonteCarloSimulators.cpp
    )
    set_target_properties(montecarlo_bindings PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_include_directories(montecarlo_bindings PRIVATE src ${MODEL_DIRS})
    target_link_libraries(montecarlo_bindings PRIVATE ${SIMULATOR_LIBS} montecarlo)
    target_compile_definitions(montecarlo_bindings PUBLIC CLUSTER_UPDATE)

    install(
        TARGETS montecarlo_bindings
        LIBRARY DESTINATION pymc
    )
endif()
