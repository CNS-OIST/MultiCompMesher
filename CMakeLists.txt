cmake_minimum_required(VERSION 3.16)
project( MultiCompMesher )

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  create_single_source_cgal_program( "src/MultiCompMesher.cpp" )
  find_package(Boost COMPONENTS program_options REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES( "MultiCompMesher" LINK_PUBLIC ${Boost_LIBRARIES} )
  else()
    message(STATUS "This program requires the Boost library, and will not be compiled.")
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
