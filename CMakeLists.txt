cmake_minimum_required( VERSION 3.0 )
project ( TAURUS VERSION 1.0 )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )
endif()

set( Boost_USE_STATIC_LIBS ON )

find_library( MATH_LIBRARIES m REQUIRED )
find_package( Threads REQUIRED )
find_package( LEMON REQUIRED )
find_package( Boost REQUIRED COMPONENTS thread system filesystem )

set( CommonLibs
  ${MATH_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${LEMON_LIBRARIES}
  ${Boost_LIBRARIES}
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories( ${LEMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} )

set( dimacs_hdr
  src/utils.h
  src/basetree.h
  src/dimacs.h
  )

set( dimacs_src
  src/dimacsmain.cpp
  src/utils.cpp
  src/basetree.cpp
  src/dimacs.cpp
  )

set( gamma_hdr
  src/utils.h
  src/basetree.h
  )

set( gamma_src
  src/utils.cpp
  src/basetree.cpp
  src/gammamain.cpp
  )

set( fpt_hdr
  src/utils.h
  src/transmissiontree.h
  src/basetree.h
  src/fptsolver.h
  )

set( fpt_src
  src/fptmain.cpp
  src/utils.cpp
  src/transmissiontree.cpp
  src/basetree.cpp
  src/fptsolver.cpp
  )

set( naive_hdr
  src/utils.h
  src/basetree.h
  src/naive.h
  )

set( naive_src
  src/naivemain.cpp
  src/utils.cpp
  src/basetree.cpp
  src/naive.cpp
  )

set( sample_naive_hdr
  src/utils.h
  src/basetree.h
  src/naive.h
  )

set( sample_naive_src
  src/samplenaivemain.cpp
  src/utils.cpp
  src/basetree.cpp
  src/naive.cpp
  )

set( sctt_hdr
  src/utils.h
  src/scttsolver.h
  )

set( sctt_src
  src/scttmain.cpp
  src/utils.cpp
  src/scttsolver.cpp
  )

set( parser_hdr
  src/utils.h
  src/basetree.h
  src/unigenparser.h
  )

set( parser_src
  src/unigenparsermain.cpp
  src/utils.cpp
  src/basetree.cpp
  src/unigenparser.cpp
  )

add_executable( dimacs ${dimacs_src} ${dimacs_hdr} )
target_link_libraries( dimacs ${CommonLibs} )

add_executable( gamma ${gamma_src} ${gamma_hdr} )
target_link_libraries( gamma ${CommonLibs} )

add_executable( fpt ${fpt_src} ${fpt_hdr} )
target_link_libraries( fpt ${CommonLibs} )

add_executable( naive ${naive_src} ${naive_hdr} )
target_link_libraries( naive ${CommonLibs} )

add_executable( naive_sample ${sample_naive_src} ${sample_naive_hdr} )
target_link_libraries( naive_sample ${CommonLibs} )

add_executable( sctt ${sctt_src} ${sctt_hdr} )
target_link_libraries( sctt ${CommonLibs} )

add_executable( unigenparser ${parser_src} ${parser_hdr} )
target_link_libraries( unigenparser ${CommonLibs} )
