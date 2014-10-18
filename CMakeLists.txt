cmake_minimum_required (VERSION 2.8)

#####
#  Credits: A significant portion of this CMakeLists.txt is adapted from cmake
#  instructions in Sailfish (Patro et al., https://github.com/kingsfordgroup/sailfish)
#####

project (MOJO)

set (MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set (MACHTYPE "x86_64")

#set version numbers
set(CPACK_PACKAGE_VERSION "0.0.2")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "2")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Chicago")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MOJO - Gene fusion discovery from RNA-seq")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

## Set the standard required compile flags
set (CMAKE_CXX_FLAGS "-fPIC -O3 -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DHAVE_CONFIG_H -Wall -std=c++11")

if (APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set (SHARED_LIB_EXTENSION "dylib")
else()
    set (SHARED_LIB_EXTENSION "so")
endif()

##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.8 or higher.")
    endif ()
    set (BOOST_CONFIGURE_TOOLSET "")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

message(STATUS "BOOST_EXTRA_FLAGS = ${BOOST_EXTRA_FLAGS}")

# Have CMake tell us what it's doing
set (CMAKE_VERBOSE_MAKEFILE true)


############
# Bwa
############
set( BWA_BINARIES
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bwa/bwa
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bwa/xa2multi.pl
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bwa/qualfa2fq.pl
)
include(ExternalProject)
message("Build system will fetch and build Bwa-0.7.10")
message("==================================================================")
ExternalProject_Add(bwa
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        URL http://dmel.uchicago.edu/~chai/MOJO/external/bwa-0.7.10.tar.bz2
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/bwa
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bwa
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/external/bwa/ && cp ${BWA_BINARIES} ${CMAKE_CURRENT_SOURCE_DIR}/external/bwa/
)

############
# bowtie2
############
include(ExternalProject)
message("Build system will fetch and build Bowtie2-2.2.3")
message("==================================================================")
set(BOWTIE2_BINARIES
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-align-l
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-align-s
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-build
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-build-l
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-build-s
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-inspect
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-inspect-l
${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2/bowtie2-inspect-s
)

ExternalProject_Add(bowtie2
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        URL http://dmel.uchicago.edu/~chai/MOJO/external/bowtie2-2.2.3-source.zip
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/bowtie2
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bowtie2
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/external/bowtie2/ && cp ${BOWTIE2_BINARIES} ${CMAKE_CURRENT_SOURCE_DIR}/external/bowtie2
)

############
# Bamtools
############
include(ExternalProject)
message("Build system will fetch and build Bamtools")
message("==================================================================")
set( CMAKE_C_FLAGS "-Wall -g -O2")
set(BAMTOOLS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/source/bamtools)
set(bamtools_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/bamtools")
set(bamtools_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${bamtools_INSTALL_DIR})

ExternalProject_Add(bamtools
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        PREFIX ${BAMTOOLS_ROOT}
        URL http://dmel.uchicago.edu/~chai/MOJO/external/BamTools.v2.3.0.tar.gz
        BINARY_DIR ${BAMTOOLS_ROOT}
        INSTALL_DIR ${bamtools_INSTALL_DIR}
        CMAKE_ARGS ${bamtools_CMAKE_ARGS}
)
set (BAM_SOURCE_DIR ${bamtools_INSTALL_DIR}/include/bamtools/ ${bamtools_INSTALL_DIR}/include/bamtools/api/ ${bamtools_INSTALL_DIR}/include/bamtools/shared/)
set (BAM_LIB_DIR ${bamtools_INSTALL_DIR}/lib/bamtools/)
set (BAM_LIBRARY ${bamtools_INSTALL_DIR}/lib/bamtools/libbamtools.${SHARED_LIB_EXTENSION})

############
# Blat
############
include(ExternalProject)
message("Build system will fetch and build Blat")
message("==================================================================")
ExternalProject_Add(Blat
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        URL http://dmel.uchicago.edu/~chai/MOJO/external/blatSrc35.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/blat
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/blat
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE} MACHTYPE=${MACHTYPE}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/external/blat/ && cp ${CMAKE_CURRENT_SOURCE_DIR}/external/source/blat/blat/blat ${CMAKE_CURRENT_SOURCE_DIR}/external/blat
)

############
# Samtools
############
include(ExternalProject)
message("Build system will fetch and build Samtools")
message("==================================================================")
set (SAMTOOLS_BINARY
${CMAKE_CURRENT_SOURCE_DIR}/external/source/samtools/samtools
)
ExternalProject_Add(samtools
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        URL http://dmel.uchicago.edu/~chai/MOJO/external/samtools-1.0.tar.bz2
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/samtools
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools/ && cp ${SAMTOOLS_BINARY} ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
)

############
# Zlib
############
include(ExternalProject)
message("Build system will fetch and build Zlib")
message("==================================================================")
set( CMAKE_C_FLAGS "-Wall -g -O2")
set(ZLIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/source/zlib)
set(zlib_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/zlib")
set(zlib_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${zlib_INSTALL_DIR})
ExternalProject_Add(ZLIB
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source
        PREFIX ${ZLIB_ROOT}
        URL http://dmel.uchicago.edu/~chai/MOJO/external/zlib-1.2.8.tar.gz
        BINARY_DIR ${ZLIB_ROOT}
        INSTALL_DIR ${zlib_INSTALL_DIR}
        CMAKE_ARGS ${zlib_CMAKE_ARGS}
)


set( BOOST_EXTRA_FLAGS "--layout=tagged" )

##
# Super-secret override
##
if ( DEFINED CUSTOM_BOOST_PATH )
        set (CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
    set (CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif ( DEFINED CUSTOM_BOOST_PATH )

##
# We want dyanmic, multithreaded boost libraries
##
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

set(Boost_ADDITIONAL_VERSIONS "1.55" "1.55.0")
find_package(Boost 1.55.0 COMPONENTS regex program_options system filesystem thread log_setup log iostreams date_time timer chrono)

# Download and install boost if not found 
include(ExternalProject)
if (NOT Boost_FOUND ) 
    message("==================================================================")
    message(STATUS
	"MOJO requires Boost 1.55.0 installed. A prior installation of boost-1.55.0 is not found (BOOST_ROOT). \n"
        "\tConfiguring to fetch and install boost \n"
    )
    set(BOOST_WILL_RECONFIGURE TRUE)
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/
        URL http://dmel.uchicago.edu/~chai/MOJO/external/boost_1_55_0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/source/boost_1_55_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55
        CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/source/boost_1_55_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} --prefix=<INSTALL_DIR>
        BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/source/boost_1_55_0/b2 ${BOOST_EXTRA_FLAGS} install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
    )

    ExternalProject_Add_Step(libboost reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} -DBOOST_WILL_RECONFUGRE=FALSE -DBOOST_RECONFIGURE=TRUE
        DEPENDEES install
    )
endif()

##
# If we're fetching boost and we need to have dummy paths for these variables
# so that CMake won't complain
##
if (BOOST_WILL_RECONFIGURE)
    message("Setting Temporary Boost paths")
    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55/include)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55/lib)
    set(Boost_FOUND TRUE)
endif()

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if (BOOST_RECONFIGURE)
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE)
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55/include)
        set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_55/lib)
        find_package(Boost 1.55.0 COMPONENTS regex program_options system filesystem thread log_setup log iostreams date_time timer chrono)
endif()


add_subdirectory ( src )

include (CPack)
