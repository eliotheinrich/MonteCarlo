cmake_minimum_required(VERSION 3.11.4)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(montecarlo
	src/MonteCarlo.h 
	src/Utility.cpp src/Utility.h
	src/DataFrame.cpp src/DataFrame.h
)
target_link_libraries(montecarlo -lpthread)
target_include_directories(montecarlo PRIVATE src)

add_subdirectory(src/Ising)
add_subdirectory(src/Clock)
add_subdirectory(src/Heisenberg)
add_subdirectory(src/XY)
add_subdirectory(src/Graph)

# Including boost packages
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Including CTPL
include_directories($ENV{CTPL_DIR})

# Including Eigen
include_directories($ENV{EIGEN_DIR})
