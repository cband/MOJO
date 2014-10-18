

message("\n")
message("Installing MOJO...")
message("==========================================================================")

FILE(MAKE_DIRECTORY ${MAIN_DIR}/lib)

message("Copying boost binaries to ${MAIN_DIR}/lib/")
FOREACH(D ${Boost_LIBRARY_DIRS}/)
        FILE( COPY ${D} DESTINATION ${MAIN_DIR}/lib/ )
ENDFOREACH()

message("Copying BamTools binaries to ${MAIN_DIR}/lib/")
FOREACH(D ${BAM_LIB_DIR}/)
        FILE( COPY ${D} DESTINATION ${MAIN_DIR}/lib/ )
ENDFOREACH()
	
message("Copying Zlib binaries to ${MAIN_DIR}/lib/")
FOREACH(D ${MAIN_DIR}/external/zlib/lib/)
        FILE( COPY ${D} DESTINATION ${MAIN_DIR}/lib/ )
ENDFOREACH()

message("Removing ${MAIN_DIR}/external/source/")
FILE(REMOVE_RECURSE ${MAIN_DIR}/external/source/)
FILE(GLOB BUILD_FILES_LIST ${MAIN_DIR}/build/)
FOREACH( D ${BUILD_FILES_LIST} )
        FILE(REMOVE_RECURSE ${D})
ENDFOREACH()
	
message("MOJO Installation complete")
if (APPLE)
        message("Please add ${MAIN_DIR}/lib to your DYLD_FALLBACK_LIBRARY_PATH")
else()
        message("Please add ${MAIN_DIR}/lib to your LD_LIBRARY_PATH")
endif()
message("==========================================================================")
