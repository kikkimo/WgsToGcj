# top-level CMake configuration file for WGSTOGCJ.ai
cmake_minimum_required(VERSION 3.12)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0079)
  cmake_policy(SET CMP0079 NEW)
endif()

project(WgsToGcj)
set(WGSTOGCJ_VERSION 1.0.0)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/sdk/lib)
set(WGSTOGCJ_WORKSPACE_DIR "$<SHELL_PATH:${CMAKE_SOURCE_DIR}/../>")
message(STATUS "WGSTOGCJ_WORKSPACE_DIR: ${WGSTOGCJ_WORKSPACE_DIR}")

option(USE_CERES "use Ceres" OFF)

# section third-party
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W1")
  endif()
  
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -DCPL_DISABLE_DLL)
  if(USE_CERES)
  	add_definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE -D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS -D_USE_CERES)
  endif()
  
  set(WGSTOGCJ_DEBUG_POSTFIX d)
  set(CDEGREL $<$<CONFIG:Debug>:debug>$<$<CONFIG:Release>:release>)
  message(STATUS "Debug_Release_Mode: ${CDEGREL}")

	set(WGSTOGCJ_BIN_ROOT "bin")
	set(WGSTOGCJ_BIN_DIR "${WGSTOGCJ_BIN_ROOT}/${CDEGREL}")
		
	if(USE_CERES)
		set(THIRD_PARTY "$ENV{WX_THIRD_PARTY}/windows")
		
		set(CERES_ROOT "${THIRD_PARTY}/ceres_2.0.0")
  	set(CERES_BIN_DIR "${CERES_ROOT}/bin")
  	set(CERES_INCLUDE_DIR "${CERES_ROOT}/include")
  	set(CERES_LIBRARY_DIR "${CERES_ROOT}/lib")
  	set(CERES_LIBS debug "${CERES_LIBRARY_DIR}/ceres-debug.lib"
                 optimized "${CERES_LIBRARY_DIR}/ceres.lib")
    
    set(GFLAGS_ROOT "${THIRD_PARTY}/gflags_2.2.2")
  	set(GFLAGS_BIN_DIR "${GFLAGS_ROOT}/bin")
  	set(GFLAGS_INCLUDE_DIR "${GFLAGS_ROOT}/include")
  	set(GFLAGS_LIBS "${GFLAGS_ROOT}/lib/gflags.lib"
                  "${GFLAGS_ROOT}/lib/gflags_nothreads.lib")

  	set(GLOG_ROOT "${THIRD_PARTY}/glog_0.4.0")
  	set(GLOG_BIN_DIR "${GLOG_ROOT}/bin")
  	set(GLOG_INCLUDE_DIR "${GLOG_ROOT}/include")
  	set(GLOG_LIBRARY "${GLOG_ROOT}/lib/glog.lib")

  	set(EIGEN_INCLUDE_DIR "${THIRD_PARTY}/eigen_3.3.90/include/eigen3")

  	set(SUITESPARSE_INCLUDE_DIR "${THIRD_PARTY}/suitesparse_metis_5.1.0/include/suitesparse")
  endif()
  
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(WGSTOGCJ_BIN_ROOT "bin")
  set(WGSTOGCJ_BIN_DIR "bin")
  set(WGSTOGCJ_DEBUG_POSTFIX "")

 	if(USE_CERES)
		set(THIRD_PARTY "$ENV{WX_THIRD_PARTY}/linux")
		
		set(CERES_ROOT "${THIRD_PARTY}/ceres_2.0.0")
  	set(CERES_BIN_DIR "${CERES_ROOT}/bin")
  	set(CERES_INCLUDE_DIR "${CERES_ROOT}/include")
  	set(CERES_LIBRARY_DIR "${CERES_ROOT}/lib")
  	set(CERES_LIBS debug "${CERES_LIBRARY_DIR}/ceres-debug.lib"
                 optimized "${CERES_LIBRARY_DIR}/ceres.lib")
    
    set(GFLAGS_ROOT "${THIRD_PARTY}/gflags_2.2.2")
  	set(GFLAGS_BIN_DIR "${GFLAGS_ROOT}/bin")
  	set(GFLAGS_INCLUDE_DIR "${GFLAGS_ROOT}/include")
  	set(GFLAGS_LIBS "${GFLAGS_ROOT}/lib/gflags.lib"
                  "${GFLAGS_ROOT}/lib/gflags_nothreads.lib")

  	set(GLOG_ROOT "${THIRD_PARTY}/glog_0.4.0")
  	set(GLOG_BIN_DIR "${GLOG_ROOT}/bin")
  	set(GLOG_INCLUDE_DIR "${GLOG_ROOT}/include")
  	set(GLOG_LIBRARY "${GLOG_ROOT}/lib/glog.lib")

  	set(EIGEN_INCLUDE_DIR "${THIRD_PARTY}/eigen_3.3.90/include/eigen3")
  endif()
                  
else()
  message(FATAL_ERROR "${CMAKE_SYSTEM_NAME}: not supported.")
endif()

message(STATUS "WGSTOGCJ_BIN_DIR: ${WGSTOGCJ_BIN_DIR}")

if(DEFINED CMAKE_DEBUG_POSTFIX)
  set(WGSTOGCJ_DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
endif()


find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(USE_CERES)
find_package(Ceres REQUIRED)
endif()



# section include
if(USE_CERES)
include_directories(${CERES_INCLUDE_DIR})
include_directories(${GFLAGS_INCLUDE_DIR})
include_directories(${GLOG_INCLUDE_DIR})
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${SUITESPARSE_INCLUDE_DIR})
endif()

include_directories("include")

# section source
add_subdirectory(src)
enable_testing()
add_subdirectory(test)

# section install
if(USE_CERES)
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  install(DIRECTORY ${CERES_BIN_DIR}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${CERES_BIN_DIR}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
  install(DIRECTORY ${GFLAGS_ROOT}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${GFLAGS_ROOT}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
  install(DIRECTORY ${GLOG_ROOT}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${GLOG_ROOT}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
else()
  # TODO
  install(DIRECTORY ${CERES_BIN_DIR}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${CERES_BIN_DIR}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
  install(DIRECTORY ${GFLAGS_ROOT}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${GFLAGS_ROOT}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
  install(DIRECTORY ${GLOG_ROOT}
          DESTINATION bin/debug
          CONFIGURATIONS Debug
          FILES_MATCHING PATTERN "*.dll" PATTERN "*.pdb")
  install(DIRECTORY ${GLOG_ROOT}
  				DESTINATION bin/release
          CONFIGURATIONS Release
          FILES_MATCHING PATTERN "*.dll")
endif()
endif()

install(DIRECTORY include DESTINATION sdk)



