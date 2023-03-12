add_library(clock
	ClockModel.h
	SquareClockModel.cpp
)
target_link_libraries(clock montecarlo)
target_include_directories(clock PRIVATE ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})