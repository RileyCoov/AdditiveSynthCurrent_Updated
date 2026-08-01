# Portable headless build for the additive resynthesis engine.
#
# The engine is pure, portable C++ (standard headers + a custom RealFFT, WAV I/O
# via <fstream>/fread, CLI: input.wav -> output.wav). It historically shipped only
# as an Xcode project; this build lets it compile and run headless on Linux/CI so
# the regression battery (Phase 0) can be reproducible off a Mac.
#
# The Xcode project remains the macOS path and is untouched.
#
#   cmake -S . -B build && cmake --build build
#   ./build/additive input.wav output.wav [args...]

cmake_minimum_required(VERSION 3.16)
project(additive_synth LANGUAGES CXX)

# gnu++17 (not strict c++17) so M_PI from <cmath> is visible on glibc, matching
# the behavior Apple clang provides by default.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB ENGINE_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/AdditiveSynthClean/*.cpp")

add_executable(additive ${ENGINE_SOURCES})
target_include_directories(additive PRIVATE "${CMAKE_SOURCE_DIR}/AdditiveSynthClean")
