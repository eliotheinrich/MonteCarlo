add_library(clock
	ClockModel.h
	SquareClockModel.h
)
target_include_directories(clock PRIVATE ${CMAKE_SOURCE_DIR}/src)
