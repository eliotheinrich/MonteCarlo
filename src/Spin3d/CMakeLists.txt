add_library(spin3d
  Spin3DModel.cpp
	TrigonalModel.cpp
  XXZHeis.cpp
  AltermagnetModel.cpp
  HelixModel.cpp
)

target_include_directories(spin3d PUBLIC ${CMAKE_SOURCE_DIR}/src)
#target_link_libraries(spin3d PUBLIC fftw3 finufft)
