cmake_minimum_required(VERSION 3.18)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 20)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_FLAGS_DEBUG "-g -funwind-tables -fasynchronous-unwind-tables")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include(FetchContent)

# Fetching externals
FetchContent_Declare(
    nlohmann
    GIT_REPOSITORY https://github.com/nlohmann/json
    GIT_TAG master
    GIT_SHALLOW TRUE
)

FetchContent_Declare(
    bspl
    GIT_REPOSITORY https://github.com/bshoshany/thread-pool
    GIT_TAG master
    GIT_SHALLOW TRUE
)

FetchContent_Declare(
    dataframe
    GIT_REPOSITORY https://github.com/eliotheinrich/dataframe
    GIT_TAG master
    GIT_SHALLOW TRUE
)

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

FetchContent_MakeAvailable(nlohmann bspl dataframe eigen)

include_directories(
    ${nlohmann_SOURCE_DIR}/include 
    ${bspl_SOURCE_DIR}/include 
    ${dataframe_SOURCE_DIR}/src
    ${eigen_SOURCE_DIR}
)

add_library(montecarlo
	src/MonteCarlo.cpp 
)
target_include_directories(montecarlo PUBLIC src)
target_link_libraries(montecarlo -lpthread)

set(ISING_DIR src/Ising)
set(CLOCK_DIR src/Clock)
set(SPIN2D_DIR src/Spin2d)
set(SPIN3D_DIR src/Spin3d)
set(GRAPH_DIR src/Graph)

set(MODEL_DIRS ${ISING_DIR} ${CLOCK_DIR} ${SPIN2D_DIR} ${SPIN3D_DIR} ${GRAPH_DIR})
foreach(DIR IN LISTS MODEL_DIRS)
	add_subdirectory(${DIR})
endforeach()

add_executable(main src/main.cpp)
set_target_properties(main PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_compile_definitions(main PUBLIC CLUSTER_UPDATE)
target_include_directories(main PRIVATE ${MODEL_DIRS})
target_link_libraries(main PRIVATE montecarlo ising spin2d spin3d clock graph)
