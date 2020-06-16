project(test)

include_directories(${CMAKE_BINARY_DIR}/test/)

option(INSTALL_TEST "install demo test" ON)

set(TESTS_BIN_DIR ${CMAKE_INSTALL_PREFIX}/${WGSTOGCJ_BIN_DIR})

if(INSTALL_TEST)
	aux_source_directory(. DIR_SRCS)
	add_executable(CovertDemo ${DIR_SRCS})
	set_target_properties(CovertDemo
    PROPERTIES
    CXX_STANDARD 17
    DEBUG_POSTFIX "${WGSTOGCJ_DEBUG_POSTFIX}")
  target_link_libraries(CovertDemo WGSTOGCJ_API)
  wgstogcj_install_target(CovertDemo)
endif()
