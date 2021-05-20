cmake_minimum_required(VERSION 3.12)
project(
        ULTRA
        VERSION 1.0.0
        HOMEPAGE_URL https://github.com/TravisWheelerLab/ULTRA
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)

set(
        LIB_CPP_FILES

        lib/json11.cpp
)
set(
        HPP_FILES

        src/FASTAReader.hpp
        src/FileReader.hpp
        src/JSONPass.hpp
        src/JSONReader.hpp
        src/JSONRepeat.hpp
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
        src/JSONPass.cpp
        src/JSONReader.cpp
        src/JSONRepeat.cpp
        src/SequenceWindow.cpp
        src/Symbol.cpp
        src/main.cpp
        src/repeat.cpp
        src/settings.cpp
        src/ultra.cpp
        src/umatrix.cpp
        src/umodel.cpp
        src/umodeldriver.cpp
)

find_package(Threads REQUIRED)

add_executable(
        ultra

        ${LIB_CPP_FILES}
        ${CPP_FILES}
)
target_link_libraries(ultra PRIVATE Threads::Threads)
install(TARGETS ultra RUNTIME)

if (NOT ULTRA_CONTAINER_VERSION)
    set(ULTRA_CONTAINER_VERSION latest)
endif ()
set(ULTRA_CONTAINER_TAG traviswheelerlab/ultra-build)

# Runs clang-format on the code
add_custom_target(
        format
        COMMAND clang-format -i ${CPP_FILES} ${HPP_FILES}
        VERBATIM
)

# Checks whether clang-format needs to be run and sets a non-zero
# exit code if so
add_custom_target(
        check-format
        COMMAND clang-format --dry-run --Werror ${CPP_FILES} ${HPP_FILES}
        VERBATIM
)

add_custom_target(
        container-build
        COMMAND docker build -f Dockerfile_build -t ${ULTRA_CONTAINER_TAG}:${ULTRA_CONTAINER_VERSION} .
        VERBATIM
)

add_custom_target(
        container-push
        COMMAND docker push ${ULTRA_CONTAINER_TAG}:${ULTRA_CONTAINER_VERSION}
        VERBATIM
)
add_dependencies(container-push container-build)

add_custom_target(
        examples
        COMMAND examples/run-all.sh
        VERBATIM
)
add_dependencies(examples ultra)

# -----------------------
# Packaging Configuration
# -----------------------

SET(ULTRA_PACKAGE_VENDOR "Travis Wheeler Lab")

# DEB (Debian / Ubuntu)

add_custom_target(
    deb-container-image
    COMMAND docker build -f Dockerfile_deb -t ultra_deb .
    VERBATIM
)

add_custom_target(
    deb-container
    COMMAND docker run --mount src=${CMAKE_SOURCE_DIR},target=/code,type=bind ultra_deb
    VERBATIM
)

add_custom_target(
    deb
    COMMAND cpack -G DEB .
    VERBATIM
)

SET(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER ${ULTRA_PACKAGE_VENDOR})
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")

# RPM

add_custom_target(
    rpm-container-image
    COMMAND docker build -f Dockerfile_rpm -t ultra_rpm .
    VERBATIM
)

add_custom_target(
    rpm-container
    COMMAND docker run --mount src=${CMAKE_SOURCE_DIR},target=/code,type=bind ultra_rpm
    VERBATIM
)

add_custom_target(
    rpm
    COMMAND cpack -G RPM .
    VERBATIM
)

SET(CPACK_RPM_FILE_NAME RPM-DEFAULT)
SET(CPACK_RPM_PACKAGE_VENDOR ${ULTRA_PACKAGE_VENDOR})

# General

SET(CPACK_PACKAGE_NAME "ultra")
SET(CPACK_PACKAGE_VENDOR "Travis Wheeler Lab")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/TravisWheelerLab/ULTRA")

SET(CPACK_GENERATOR "DEB")
INCLUDE(CPack)

