add_library(xymodel
	XYModel.cpp XYModel.h
	SquareXYModel.cpp SquareXYModel.h
	TrigonalXYModel.cpp TrigonalXYModel.h
)
target_link_libraries(xymodel montecarlo)
target_include_directories(xymodel PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})