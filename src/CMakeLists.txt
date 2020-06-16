project(src)

include_directories(${CMAKE_BINARY_DIR}/src/)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

add_definitions(-D_WGSTOGCJ_API)
aux_source_directory(. DIR_SRCS)
add_library(WGSTOGCJ_API SHARED ${DIR_SRCS})
set_target_properties(WGSTOGCJ_API
    PROPERTIES
    CXX_STANDARD 17
    DEBUG_POSTFIX "${WGSTOGCJ_DEBUG_POSTFIX}")
    
if(USE_CERES)
	 target_link_libraries(WGSTOGCJ_API ${CERES_LIBS})
	 target_link_libraries(WGSTOGCJ_API ${GFLAGS_LIBS})
	 target_link_libraries(WGSTOGCJ_API ${GLOG_LIBRARY})
endif()
    
function(wgstogcj_install_target name)
  install(TARGETS ${name} ${ARGN}
    LIBRARY DESTINATION sdk/lib
    ARCHIVE DESTINATION sdk/lib
    RUNTIME DESTINATION ${WGSTOGCJ_BIN_DIR})
  if(WIN32 AND MSVC AND (NOT MSVC_VERSION LESS 1600))
    install(FILES $<TARGET_PDB_FILE:${name}> DESTINATION ${WGSTOGCJ_BIN_DIR} OPTIONAL)
  endif()
endfunction()

wgstogcj_install_target(WGSTOGCJ_API)



