cmake_minimum_required(VERSION 3.12)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 20)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

set(DF_DIR $ENV{DF_DIR})
set(NLOHMANN_DIR $ENV{NLOHMANN_DIR})
set(CTPL_DIR $ENV{CTPL_DIR})
set(EIGEN_DIR $ENV{EIGEN_DIR})

set(DF_DEPENDENCIES ${DF_DIR} ${NLOHMANN_DIR} ${CTPL_DIR})

add_library(montecarlo
	src/MonteCarlo.cpp 
)
target_include_directories(montecarlo PUBLIC src ${DF_DEPENDENCIES})
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

# Including boost packages
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(main src/main.cpp)
set_target_properties(main PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_compile_definitions(main PUBLIC CLUSTER_UPDATE)
target_include_directories(main PRIVATE 
	${DF_DEPENDENCIES} 
	${MODEL_DIRS}
)
target_link_libraries(main PRIVATE montecarlo ising spin2d spin3d clock graph)

add_executable(metropolis src/main.cpp)
set_target_properties(metropolis PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_include_directories(metropolis PRIVATE
	${DF_DEPENDENCIES}
	${MODEL_DIRS}
)
target_link_libraries(metropolis PRIVATE montecarlo ising_metropolis spin2d_metropolis spin3d_metropolis clock graph)

