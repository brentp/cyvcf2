cmake_minimum_required(VERSION 3.20...3.30)

project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

if(SKBUILD_STATE STREQUAL "sdist")
  find_program(MAKE_EXECUTABLE NAMES make gmake mingw32-make REQUIRED)
  execute_process(
    COMMAND "${MAKE_EXECUTABLE}" htscodecs/htscodecs/version.h
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/htslib"
    COMMAND_ERROR_IS_FATAL ANY
  )
  return()
endif()

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

set(CYVCF2_HTSLIB_CONFIGURE_OPTIONS "$ENV{CYVCF2_HTSLIB_CONFIGURE_OPTIONS}" CACHE STRING "Extra options passed to htslib configure")
set(CYVCF2_HTSLIB_MODE "$ENV{CYVCF2_HTSLIB_MODE}" CACHE STRING "How to link htslib: BUILTIN or EXTERNAL")
set(CYVCF2_CYTHONIZE "$ENV{CYTHONIZE}" CACHE STRING "Whether to regenerate cyvcf2.c from Cython: AUTO, ON, or OFF")

if(CYVCF2_HTSLIB_MODE STREQUAL "")
  if(WIN32)
    set(CYVCF2_HTSLIB_MODE "EXTERNAL")
  else()
    set(CYVCF2_HTSLIB_MODE "BUILTIN")
  endif()
endif()

string(TOUPPER "${CYVCF2_HTSLIB_MODE}" CYVCF2_HTSLIB_MODE)
set_property(CACHE CYVCF2_HTSLIB_MODE PROPERTY STRINGS BUILTIN EXTERNAL)

if(NOT CYVCF2_HTSLIB_MODE STREQUAL "BUILTIN" AND NOT CYVCF2_HTSLIB_MODE STREQUAL "EXTERNAL")
  message(FATAL_ERROR "CYVCF2_HTSLIB_MODE must be BUILTIN or EXTERNAL, got '${CYVCF2_HTSLIB_MODE}'")
endif()

set(CYVCF2_SOURCE_C "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2/cyvcf2.c")
set(CYVCF2_GENERATED_C "${CMAKE_CURRENT_BINARY_DIR}/generated/cyvcf2.c")

function(cyvcf2_boolish INPUT OUTPUT)
  string(TOUPPER "${INPUT}" NORMALIZED_INPUT)
  if(NORMALIZED_INPUT STREQUAL "" OR NORMALIZED_INPUT STREQUAL "AUTO")
    set(${OUTPUT} "AUTO" PARENT_SCOPE)
  elseif(NORMALIZED_INPUT MATCHES "^(1|ON|TRUE|YES)$")
    set(${OUTPUT} "ON" PARENT_SCOPE)
  elseif(NORMALIZED_INPUT MATCHES "^(0|OFF|FALSE|NO)$")
    set(${OUTPUT} "OFF" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Expected a boolean-like value or AUTO, got '${INPUT}'")
  endif()
endfunction()

cyvcf2_boolish("${CYVCF2_CYTHONIZE}" CYVCF2_CYTHONIZE_NORMALIZED)

if(CYVCF2_CYTHONIZE_NORMALIZED STREQUAL "AUTO")
  if(EXISTS "${CYVCF2_SOURCE_C}")
    set(CYVCF2_CYTHONIZE_NORMALIZED "OFF")
  else()
    set(CYVCF2_CYTHONIZE_NORMALIZED "ON")
  endif()
endif()

if(CYVCF2_CYTHONIZE_NORMALIZED)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
  add_custom_command(
    OUTPUT "${CYVCF2_GENERATED_C}"
    COMMAND Python::Interpreter -m cython
            "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2/cyvcf2.pyx"
            --output-file "${CYVCF2_GENERATED_C}"
    DEPENDS
      "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2/cyvcf2.pyx"
      "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2/cyvcf2.pxd"
    COMMENT "Generating cyvcf2.c from Cython"
    VERBATIM
  )
  set(CYVCF2_EXTENSION_SOURCE "${CYVCF2_GENERATED_C}")
else()
  if(NOT EXISTS "${CYVCF2_SOURCE_C}")
    message(FATAL_ERROR "CYTHONIZE is OFF, but ${CYVCF2_SOURCE_C} does not exist")
  endif()
  set(CYVCF2_EXTENSION_SOURCE "${CYVCF2_SOURCE_C}")
endif()

function(cyvcf2_split_configure_options OUTPUT)
  if(CYVCF2_HTSLIB_CONFIGURE_OPTIONS STREQUAL "")
    set(${OUTPUT} "" PARENT_SCOPE)
  else()
    separate_arguments(_cyvcf2_configure_options UNIX_COMMAND "${CYVCF2_HTSLIB_CONFIGURE_OPTIONS}")
    set(${OUTPUT} "${_cyvcf2_configure_options}" PARENT_SCOPE)
  endif()
endfunction()

function(cyvcf2_build_htslib STATIC_MODE)
  find_program(MAKE_EXECUTABLE NAMES make gmake mingw32-make REQUIRED)

  set(HTSLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
  cyvcf2_split_configure_options(HTSLIB_CONFIGURE_OPTION_LIST)

  if(NOT EXISTS "${HTSLIB_SOURCE_DIR}/config.status")
    execute_process(
      COMMAND autoreconf -i
      WORKING_DIRECTORY "${HTSLIB_SOURCE_DIR}"
      COMMAND_ERROR_IS_FATAL ANY
    )

    set(HTSLIB_CONFIGURE_COMMAND ./configure)
    if(STATIC_MODE)
      list(APPEND HTSLIB_CONFIGURE_COMMAND "CFLAGS=-fPIC")
    endif()
    list(APPEND HTSLIB_CONFIGURE_COMMAND ${HTSLIB_CONFIGURE_OPTION_LIST})

    execute_process(
      COMMAND ${HTSLIB_CONFIGURE_COMMAND}
      WORKING_DIRECTORY "${HTSLIB_SOURCE_DIR}"
      COMMAND_ERROR_IS_FATAL ANY
    )
  endif()

  execute_process(
    COMMAND "${MAKE_EXECUTABLE}" htscodecs/htscodecs/version.h
    WORKING_DIRECTORY "${HTSLIB_SOURCE_DIR}"
    COMMAND_ERROR_IS_FATAL ANY
  )

  if(STATIC_MODE)
    set(HTSLIB_BUILD_TARGET libhts.a)
  elseif(APPLE)
    set(HTSLIB_BUILD_TARGET libhts.dylib)
  elseif(WIN32)
    set(HTSLIB_BUILD_TARGET hts-3.dll)
  else()
    set(HTSLIB_BUILD_TARGET libhts.so)
  endif()

  execute_process(
    COMMAND "${MAKE_EXECUTABLE}" "${HTSLIB_BUILD_TARGET}"
    WORKING_DIRECTORY "${HTSLIB_SOURCE_DIR}"
    COMMAND_ERROR_IS_FATAL ANY
  )
endfunction()

function(cyvcf2_read_htslib_static_libs OUTPUT)
  set(HTSLIB_CONFIG_STATUS "${CMAKE_CURRENT_SOURCE_DIR}/htslib/config.status")
  if(NOT EXISTS "${HTSLIB_CONFIG_STATUS}")
    set(${OUTPUT} "" PARENT_SCOPE)
    return()
  endif()

  file(READ "${HTSLIB_CONFIG_STATUS}" HTSLIB_CONFIG_STATUS_CONTENTS)
  string(REGEX MATCH "S\\[\"static_LIBS\"\\]=\"([^\"]*)\"" _static_libs_match "${HTSLIB_CONFIG_STATUS_CONTENTS}")
  if(_static_libs_match)
    separate_arguments(HTSLIB_STATIC_LIBS UNIX_COMMAND "${CMAKE_MATCH_1}")
    set(${OUTPUT} "${HTSLIB_STATIC_LIBS}" PARENT_SCOPE)
  else()
    set(${OUTPUT} "" PARENT_SCOPE)
  endif()
endfunction()

set(CYVCF2_HTSLIB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/htslib" "${CMAKE_CURRENT_SOURCE_DIR}/htslib/htslib")
set(CYVCF2_HTSLIB_LIBRARIES "")
set(CYVCF2_HTSLIB_TARGET "")

if(CYVCF2_HTSLIB_MODE STREQUAL "BUILTIN")
  if(WIN32)
    message(FATAL_ERROR "CYVCF2_HTSLIB_MODE=BUILTIN is not supported on Windows")
  endif()

  message(STATUS "cyvcf2: building bundled htslib")
  cyvcf2_build_htslib(TRUE)
  cyvcf2_read_htslib_static_libs(CYVCF2_HTSLIB_STATIC_LIBS)

  set(CYVCF2_HTSLIB_STATIC_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a")
  if(NOT EXISTS "${CYVCF2_HTSLIB_STATIC_LIBRARY}")
    message(FATAL_ERROR "Expected ${CYVCF2_HTSLIB_STATIC_LIBRARY} after building htslib")
  endif()

  add_library(cyvcf2_htslib STATIC IMPORTED GLOBAL)
  set_target_properties(cyvcf2_htslib PROPERTIES
    IMPORTED_LOCATION "${CYVCF2_HTSLIB_STATIC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CYVCF2_HTSLIB_INCLUDE_DIRS}"
  )
  set(CYVCF2_HTSLIB_TARGET cyvcf2_htslib)
  list(APPEND CYVCF2_HTSLIB_LIBRARIES ${CYVCF2_HTSLIB_STATIC_LIBS})
else()
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    if(WIN32)
      pkg_check_modules(HTSLIB_PC QUIET htslib)
    else()
      pkg_check_modules(HTSLIB_PC QUIET IMPORTED_TARGET htslib)
    endif()
  endif()

  if(HTSLIB_PC_FOUND)
    if(WIN32)
      find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
      if(NOT HTSLIB_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find external htslib headers")
      endif()
      list(APPEND CYVCF2_HTSLIB_INCLUDE_DIRS "${HTSLIB_INCLUDE_DIR}")
      list(APPEND CYVCF2_HTSLIB_LIBRARIES ${HTSLIB_PC_LIBRARIES})
    else()
      list(APPEND CYVCF2_HTSLIB_LIBRARIES PkgConfig::HTSLIB_PC)
    endif()
  else()
    find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
    find_library(HTSLIB_LIBRARY NAMES hts hts-3 libhts)
  endif()

  if(NOT HTSLIB_PC_FOUND AND (NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY))
    if(WIN32)
      message(FATAL_ERROR "Could not find external htslib; install htslib first")
    endif()

    message(STATUS "cyvcf2: external htslib not found; building bundled shared htslib")
    cyvcf2_build_htslib(FALSE)
    set(HTSLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
    find_library(HTSLIB_LIBRARY NAMES hts PATHS "${CMAKE_CURRENT_SOURCE_DIR}/htslib" NO_DEFAULT_PATH)
  endif()

  if(NOT HTSLIB_PC_FOUND AND NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR "Could not locate htslib library")
  endif()

  if(NOT HTSLIB_PC_FOUND)
    list(APPEND CYVCF2_HTSLIB_INCLUDE_DIRS "${HTSLIB_INCLUDE_DIR}")
    list(APPEND CYVCF2_HTSLIB_LIBRARIES "${HTSLIB_LIBRARY}")
  endif()
endif()

Python_add_library(cyvcf2_extension MODULE WITH_SOABI
  "${CYVCF2_EXTENSION_SOURCE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2/helpers.c"
)

set_target_properties(cyvcf2_extension PROPERTIES
  OUTPUT_NAME "cyvcf2"
)

target_include_directories(cyvcf2_extension PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/cyvcf2"
  ${CYVCF2_HTSLIB_INCLUDE_DIRS}
  ${Python_NumPy_INCLUDE_DIRS}
)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(cyvcf2_extension PRIVATE
    -Wno-sign-compare
    -Wno-unused-function
    -Wno-strict-prototypes
    -Wno-unused-result
  )
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  target_compile_options(cyvcf2_extension PRIVATE -Wno-discarded-qualifiers)
endif()

if(CYVCF2_HTSLIB_TARGET)
  target_link_libraries(cyvcf2_extension PRIVATE ${CYVCF2_HTSLIB_TARGET})
endif()
target_link_libraries(cyvcf2_extension PRIVATE ${CYVCF2_HTSLIB_LIBRARIES})

install(TARGETS cyvcf2_extension
  LIBRARY DESTINATION cyvcf2
  RUNTIME DESTINATION cyvcf2
)
