
set ( FILTER_MAIN_SRCS FilterJunctAlignOutput.cpp Utils.cpp )
set ( STREAM_MAIN_SRCS StreamNthFastqSplit.cpp )

set ( MOJO_MAIN_SRCS
MOJO.cpp
Config.cpp
DiscordantReadFinder.cpp
DiscordantClusterFinder.cpp
FastqParser.cpp
FusionCompiler.cpp
FusionJunction.cpp
FusionQuant.cpp
GeneModel.cpp
GeneModelObjs.cpp
JunctionAligner.cpp
JunctionFilter.cpp
Logger.cpp
MOJO.cpp
PSLParser.cpp
Read.cpp
Utils.cpp
)

include_directories(
${MAIN_DIR}/include
${Boost_INCLUDE_DIRS}
${BAM_SOURCE_DIR}
)

link_directories(
${MAIN_DIR}/lib
${Boost_LIBRARY_DIRS}
${BAM_LIB_DIR}
)

message ("lib: ${BAM_LIB_DIR}")
message ("src: ${BAM_SOURCE_DIR}")
message ("libr: ${BAM_LIBRARY}")

add_executable( StreamNthFastqSplit StreamNthFastqSplit.cpp )
add_executable( MOJO ${MOJO_MAIN_SRCS} )
add_executable( FilterJunctAlignOutput ${FILTER_MAIN_SRCS} )

target_link_libraries( FilterJunctAlignOutput ${Boost_LIBRARIES} )
target_link_libraries( MOJO ${Boost_LIBRARIES} ${BAM_LIBRARY} )

FILE(MAKE_DIRECTORY ${MAIN_DIR}/bin)

add_custom_command(
        TARGET MOJO POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/MOJO ${MAIN_DIR}/bin/
)

add_custom_command(
        TARGET StreamNthFastqSplit POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/StreamNthFastqSplit ${MAIN_DIR}/bin/
)

add_custom_command(
        TARGET FilterJunctAlignOutput POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/FilterJunctAlignOutput ${MAIN_DIR}/bin/
)


#http://stackoverflow.com/questions/22006908/cmake-how-to-execute-a-command-before-make-install
#ADD_CUSTOM_TARGET(distclean COMMAND ${CMAKE_COMMAND} -P ${MAIN_DIR}/cmake/CleanUp.cmake)
#set(CLEANUP_SCRIPT /lustre/beagle/cbandlam/MOJO-P/cmake/CleanUp.cmake)
install(
    CODE
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
                            -DBoost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}
                            -DMAIN_DIR=${MAIN_DIR}
                            -DBAM_LIB_DIR=${BAM_LIB_DIR}
                            -P \"${MAIN_DIR}/cmake/PostInstall.cmake\")
    "
)
