cmake_minimum_required(VERSION 3.12)
project(
        ULTRA
        HOMEPAGE_URL https://github.com/TravisWheelerLab/ULTRA
)

# By default we do NOT build the standalone ULTRA library
option(BUILD_ULTRA_LIB "Also build ULTRA as a standalone library" OFF)

# Default to Release if no build type is set
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "Setting CMAKE_BUILD_TYPE=Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Add a debug flag when not in Release
if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DDEBUG_PRAGMA=1)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header files (for installation / PUBLIC_HEADER)
set(HPP_FILES
        src/FASTAReader.hpp
        src/FileReader.hpp
        src/SequenceWindow.hpp
        src/Symbol.hpp
        src/repeat.hpp
        src/ultra.hpp
        src/umatrix.hpp
        src/umodel.hpp
        src/umodeldriver.hpp
        src/BEDFileWriter.hpp
        src/TabFileWriter.hpp
        src/RepeatFileWriter.hpp
        src/JSONFileWriter.hpp
        src/RepeatSplitter.hpp
        src/cli.hpp
        src/mask.hpp
        )

# All source files
set(CPP_FILES
        src/FASTAReader.cpp
        src/FileReader.cpp
        src/SequenceWindow.cpp
        src/Symbol.cpp
        src/main.cpp
        src/repeat.cpp
        src/ultra.cpp
        src/umatrix.cpp
        src/umodel.cpp
        src/umodeldriver.cpp
        src/BEDFileWriter.cpp
        src/TabFileWriter.cpp
        src/JSONFileWriter.cpp
        src/RepeatSplitter.cpp
        src/mask.cpp
        src/cli.cpp
        )

find_package(Threads REQUIRED)

# LIB_SOURCES = everything except main.cpp
set(LIB_SOURCES ${CPP_FILES})
list(REMOVE_ITEM LIB_SOURCES src/main.cpp)

# 1) Build static library only if requested
if (BUILD_ULTRA_LIB)
    add_library(ultra_core STATIC ${LIB_SOURCES})
    target_include_directories(ultra_core
            PUBLIC ${CMAKE_SOURCE_DIR}/src
            )
    target_link_libraries(ultra_core
            PUBLIC Threads::Threads
            )

    install(
            TARGETS ultra_core
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION include
    )
endif()

# 2) Build the CLI executable
if (BUILD_ULTRA_LIB)
    # Link only main.cpp against ultra_core
    add_executable(ultra src/main.cpp)
    target_link_libraries(ultra PRIVATE ultra_core)
else()
    # Monolithic: compile everything into the exe
    add_executable(ultra ${CPP_FILES})
    target_link_libraries(ultra PRIVATE Threads::Threads)
endif()

target_include_directories(ultra PRIVATE ${CMAKE_SOURCE_DIR}/src)

install(TARGETS ultra RUNTIME DESTINATION bin)

# Optional examples target
add_custom_target(
        examples
        COMMAND examples/run-all.sh
        VERBATIM
)
add_dependencies(examples ultra)
