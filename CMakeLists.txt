cmake_minimum_required(VERSION 3.0)

project(cquadpack C)

include_directories(include)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(source_files
    src/dqag.c
    src/dqage.c
    src/dqagi.c
    src/dqagp.c
    src/dqags.c
    src/dqawc.c
    src/dqawce.c
    src/dqawf.c
    src/dqawfe.c
    src/dqawo.c
    src/dqaws.c
    src/dqawse.c
    src/dqc25c.c
    src/dqc25o.c
    src/dqc25s.c
    src/dqcheb.c
    src/dqext.c
    src/dqfour.c
    src/dqk15.c
    src/dqk15i.c
    src/dqk15w.c
    src/dqk21.c
    src/dqk31.c
    src/dqk41.c
    src/dqk51.c
    src/dqk61.c
    src/dqmomo.c
    src/dqng.c
    src/dqsort.c
    src/dqwgt.c
    src/logdiff.c
    src/logsum.c
)
set(header_files
    include/cquadpack.h
)
set(BUILD_SHARED_LIBS ON)
add_library(cquadpack ${source_files} ${header_files})
message(STATUS "Will compile the library in shared mode")
set_property(TARGET cquadpack PROPERTY POSITION_INDEPENDENT_CODE 1)

# Add export header
include(GenerateExportHeader)
generate_export_header(cquadpack)
include_directories(${PROJECT_BINARY_DIR})

# Add the math libraries (gcc requires this)
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(cquadpack m)
endif()

set_target_properties(cquadpack PROPERTIES SUFFIX ".so")
install(TARGETS cquadpack DESTINATION ${CMAKE_SOURCE_DIR}/NumbaQuadpack/)

