add_library(clock
	ClockModel.h
	SquareClockModel.cpp SquareClockModel.h
)
target_link_libraries(clock montecarlo)
target_include_directories(clock PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})