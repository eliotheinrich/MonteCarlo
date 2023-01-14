cmake_minimum_required(VERSION 3.25.1)
project(MonteCarlo)

set(CMAKE_CXX_STANDARD 17)

set(PROJECT_DIR ${PROJECT_SOURCE_DIR})

set(HEADER_FILES ${PROJECT_DIR}/DataFrame.h ${PROJECT_DIR}/MonteCarlo.h)
set(SOURCE ${PROJECT_DIR}/DataFrame.cpp ${PROJECT_DIR}/MonteCarlo.cpp ${PROJECT_DIR}/ctpl.h)


add_executable(test test.cpp ${SOURCE})
target_include_directories(test PRIVATE ${PROJECT_DIR})


# Including boost packages
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
target_link_libraries(test LINK_PUBLIC ${Boost_LIBRARIES})