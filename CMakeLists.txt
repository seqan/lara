# ===========================================================================
#               LaRA -- Lagrangian Relaxed structural Alignment
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

# ----------------------------------------------------------------------------
# Compiler detection
# ----------------------------------------------------------------------------
message ("${ColourBold}Compiler detection${ColourReset}")

# change version after every release
project(lara
        VERSION 2.0.1
        LANGUAGES CXX)
set(PROJECT_DESCRIPTION "LaRA: Lagrangian Relaxed structural Alignment")
set (LARA_SOURCE_FILES src/lara.cpp)
add_executable (lara ${LARA_SOURCE_FILES})

# Make "Release" the default cmake build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()
set(CMAKE_CXX_STANDARD 14)

# Compiler specifics: Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ftemplate-depth-1024")

    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla-extension")
endif ()

# Compiler specifics: GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla")

    # parallelize parts of build even for one translation unit
    if (NOT DEFINED LARA_COMPILE_THREADS)
        include(ProcessorCount)
        ProcessorCount(LARA_COMPILE_THREADS)
    endif ()
    if (LARA_COMPILE_THREADS GREATER 1)
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto=${LARA_COMPILE_THREADS}")
    endif()

    # strip binaries to make them smaller
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif ()

# Compiler specifics: All but Intel
if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # -Wextra not met by current SeqAn on intel
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wextra")
endif ()

# ----------------------------------------------------------------------------
# Dependency detection
# ----------------------------------------------------------------------------
message ("\n${ColourBold}Dependency detection${ColourReset}")

# SeqAn
set (MINIMUM_SEQAN_VERSION "2.4.0")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
endif ()
find_package(SeqAn REQUIRED CONFIG)

# LEMON
set (LEMON_FOUND "FALSE")
find_path(LEMON_PATH NAMES lemon/smart_graph.h HINTS ${CMAKE_CURRENT_SOURCE_DIR}/include/lemon/include)
if (LEMON_PATH)
    set (CMAKE_INCLUDE_PATH ${LEMON_PATH} ${CMAKE_INCLUDE_PATH})
    add_definitions(-DLEMON_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLEMON_ONLY_TEMPLATES")
    set (LEMON_FOUND "TRUE")
else()
    message (STATUS "Cannot find header files of the LEMON library, which is required to compute matchings.")
endif ()

# Include MPFR in case ViennaRNA is compiled with it.
find_library (MPFR mpfr)
if (MPFR)
    target_link_libraries (lara PUBLIC ${MPFR})
endif ()

#ViennaRNA
set (VIENNA_FOUND "FALSE")
find_library (VIENNA_RNA_LIB libRNA.a)
find_path (VIENNA_RNA_PATH NAMES ViennaRNA/part_func.h)
if (VIENNA_RNA_LIB AND VIENNA_RNA_PATH)
    add_definitions (-DVIENNA_RNA_FOUND)
    set (VIENNA_FOUND "TRUE")
    target_link_libraries (lara PUBLIC ${VIENNA_RNA_LIB})
    target_include_directories (lara SYSTEM PUBLIC ${VIENNA_RNA_PATH})
endif ()

# others
find_package(OpenMP QUIET)
find_package(ZLIB   QUIET)
find_package(BZip2  QUIET)

message(STATUS "These dependencies where found:")
message(   "     SEQAN       ${SEQAN_FOUND}   \t${SEQAN_VERSION_STRING}")
message(   "     ZLIB        ${ZLIB_FOUND}   \t${ZLIB_VERSION_STRING}")
message(   "     BZIP2       ${BZIP2_FOUND}   \t${BZIP2_VERSION_STRING}")
message(   "     OPENMP      ${OPENMP_FOUND}   \t${OpenMP_CXX_FLAGS}")
message(   "     ViennaRNA   ${VIENNA_FOUND}   \t${VIENNA_RNA_PATH}")
message(   "     Lemon       ${LEMON_FOUND}   \t${LEMON_PATH}")

# Warn if OpenMP was not found.
if (OPENMP_FOUND)
    add_definitions(-DWITH_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else ()
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. LaRA will be built without multi-threading! "
    "This is probably not what you want! Use GCC >= 4.9.1, Clang >= 3.8.0 or ICC >= 16.0.2\nWARNING WARNING WARNING")
endif (OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (WARNING "WARNING: Zlib not found. Building LaRA without support for gzipped input and output.")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (WARNING "WARNING: BZip2 not found. Building LaRA without support for bzipped input and output.")
endif (NOT BZIP2_FOUND)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
        return ()
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        message (FATAL_ERROR "Your Clang version is too old. Please upgrade to 3.8.0 or use GCC.")
        return()
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
        message(WARNING "Your Intel Compiler version is too old. Please upgrade to 17.0.0 or newer!")
    endif ()
else ()
    message(WARNING "Unknown compiler, you are on your own!")
endif ()

# SeqAn version check
if (SEQAN_VERSION_STRING VERSION_LESS "${MINIMUM_SEQAN_VERSION}")
    message (FATAL_ERROR "The minimum SeqAn version required is ${MINIMUM_SEQAN_VERSION}!")
    return ()
endif ()

# deactivate the version check on old seqan releases
if (SEQAN_VERSION_STRING VERSION_LESS "2.3.2")
    add_definitions (-DSEQAN_DISABLE_VERSION_CHECK="YES")
endif ()

set (SEQAN_APP_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set (SEQAN_APP_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set (SEQAN_APP_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set (SEQAN_APP_VERSION "${PROJECT_VERSION}")

target_link_libraries (lara PUBLIC ${SEQAN_LIBRARIES})
message(STATUS "The requirements where met.")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Enable global exception handler for all seqan apps.
set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1)

# Add include directories.
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS} ${LEMON_PATH})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Set the right output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic")

# Warn if cmake build type is not "Release"
if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

# message ("\n${ColourBold}Setting up unit tests${ColourReset}")
# add_subdirectory(tests)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# Information
set (CPACK_PACKAGE_NAME "lara")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "LaRA: Lagrangian Relaxed structural Alignment")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VENDOR "Jörg Winkler <j.winkler@fu-berlin.de>")
set (CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
set (CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")

# Package format(s)
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CPACK_GENERATOR "ZIP;NSIS")
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CPACK_GENERATOR "ZIP;DragNDrop")
elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
    set(CPACK_GENERATOR "TBZ2")
else()
    set(CPACK_GENERATOR "TXZ")
endif ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
endif ()

# Package architecture
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
    set(CMAKE_SYSTEM_PROCESSOR "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
endif ()
if (CMAKE_CXX_FLAGS MATCHES "avx2")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_avx2")
elseif (CMAKE_CXX_FLAGS MATCHES "sse4")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_sse4")
endif()

# Include architecture in package name
if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

include (CPack)
