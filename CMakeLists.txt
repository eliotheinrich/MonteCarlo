cmake_minimum_required(VERSION 3.11.4)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -DCLUSTER_UPDATE")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DCLUSTER_UPDATE")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

set(DF_DIR $ENV{DF_DIR})
set(NLOHMANN_DIR $ENV{NLOHMANN_DIR})
set(CTPL_DIR $ENV{CTPL_DIR})
set(EIGEN_DIR $ENV{EIGEN_DIR})

set(DF_DEPENDENCIES ${DF_DIR} ${NLOHMANN_DIR} ${CTPL_DIR})
message(${DF_DEPENDENCIES})

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
target_include_directories(main PRIVATE 
	${DF_DEPENDENCIES} 
	${MODEL_DIRS}
)
target_link_libraries(main PRIVATE montecarlo ising spin2d spin3d clock graph)