cmake_minimum_required(VERSION 3.15)

project(qsystem)

if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
                  "${CMAKE_BINARY_DIR}/conan.cmake")
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)

conan_cmake_run(REQUIRES boost/1.72.0
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/armadillo-code/include)

file(GLOB_RECURSE SRC src/*.cpp src/*.c)

add_library(qsystem SHARED
            ${SRC})
            
SET_TARGET_PROPERTIES(qsystem PROPERTIES 
                      CXX_STANDARD 17
                      PREFIX ""
                      OUTPUT_NAME "_qsystem")

TARGET_LINK_LIBRARIES(qsystem ${PYTHON_LIBRARIES} CONAN_PKG::boost)

install(TARGETS qsystem DESTINATION qsystem/)
