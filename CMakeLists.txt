cmake_minimum_required(VERSION 3.15)

project(MonteCarlo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Optional: specify output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist)

# Include source directories
include_directories(
    src
    /opt/homebrew/anaconda3/envs/plots/include/eigen3
)

# Add sources
add_executable(main
    src/main.cpp
    src/Spin2d/Spin2dModel.cpp
    src/Spin2d/SquareXYModel.cpp
    src/Spin3d/Spin3dModel.cpp
    src/Ising/IsingModel.cpp
    src/Ising/SquareIsingModel.cpp
    src/MonteCarlo.cpp
    src/SimulatorDisplay/Animator.cpp
    src/SimulatorDisplay/Shader.cpp
    src/glad.cpp
)

# Compiler options
target_compile_options(main PRIVATE
    -Wall
    -Wno-unused-function
    -Wno-deprecated-declarations
)

if (DEFINED EMSCRIPTEN)
  # Force Emscripten to output HTML instead of JS
  set(CMAKE_EXECUTABLE_SUFFIX ".html")

  # Emscripten linker options
  target_link_options(main PRIVATE
    "-sUSE_GLFW=3"
    "-sFULL_ES2=1"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sEXPORTED_FUNCTIONS=['_main']"
    "-sEXPORTED_RUNTIME_METHODS=['cwrap','ccall']"
    "-sMODULARIZE=0"
    "-sEXPORT_ES6=0"
    )
  message(STATUS "Executable output: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/main.html")
else()
  target_include_directories(main PRIVATE "/opt/homebrew/include/")
  target_link_directories(main PRIVATE "/opt/homebrew/lib/")
  target_link_libraries(main PRIVATE glfw)
endif()
