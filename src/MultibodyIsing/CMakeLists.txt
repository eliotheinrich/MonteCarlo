add_library(multibody_ising
    MultibodyIsingModel.cpp
    LDPCIsingModel.cpp
)
target_include_directories(multibody_ising PUBLIC ${CMAKE_SOURCE_DIR}/src)
