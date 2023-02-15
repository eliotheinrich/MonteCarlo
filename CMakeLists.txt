cmake_minimum_required(VERSION 3.11.4)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(montecarlo
	src/ctpl.h
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

add_executable(test src/test.cpp)
target_include_directories(test PRIVATE src src/Ising)

# Including boost packages
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
target_link_libraries(test LINK_PUBLIC 
	montecarlo
	ising
	${Boost_LIBRARIES}
)

# Including Eigen
include_directories($ENV{EIGEN_DIR})
