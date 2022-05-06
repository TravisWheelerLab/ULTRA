cmake_minimum_required(VERSION 3.12)
project(
        ULTRA
        HOMEPAGE_URL https://github.com/TravisWheelerLab/ULTRA
)

set(CMAKE_CXX_STANDARD 17)
set(
        LIB_CPP_FILES

        lib/json11.cpp
)
set(
        HPP_FILES

        src/FASTAReader.hpp
        src/FileReader.hpp
        src/SequenceWindow.hpp
        src/Symbol.hpp
        src/repeat.hpp
        src/settings.hpp
        src/ultra.hpp
        src/umatrix.hpp
        src/umodel.hpp
        src/umodeldriver.hpp

)
set(
        CPP_FILES

        src/FASTAReader.cpp
        src/FileReader.cpp
        src/SequenceWindow.cpp
        src/Symbol.cpp
        src/main.cpp
        src/repeat.cpp
        src/settings.cpp
        src/ultra.cpp
        src/umatrix.cpp
        src/umodel.cpp
        src/umodeldriver.cpp
        src/BEDFileWriter.cpp src/BEDFileWriter.hpp src/RepeatFileWriter.hpp src/JSONFileWriter.cpp src/JSONFileWriter.hpp src/RepeatSplitter.cpp src/RepeatSplitter.h)

find_package(Threads REQUIRED)

add_executable(
        ultra

        ${LIB_CPP_FILES}
        ${CPP_FILES}
)
target_link_libraries(ultra PRIVATE Threads::Threads)
install(TARGETS ultra RUNTIME)

add_custom_target(
        examples
        COMMAND examples/run-all.sh
        VERBATIM
)
add_dependencies(examples ultra)
