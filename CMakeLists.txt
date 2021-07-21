cmake_minimum_required(VERSION 3.16)
project(MultiCompMesher)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
find_package(
  CGAL
  COMPONENTS ImageIO
  REQUIRED)

option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if(CGAL_FOUND)
  if(CGAL_MAJOR_VERSION LESS 5)
    message(SEND_ERROR "CGAL ${CGAL_VERSION} was found but version 5 or higher is required.")
  endif()
  include(${CGAL_USE_FILE})
  create_single_source_cgal_program("src/MultiCompMesher.cpp")
  find_package(
    Boost
    COMPONENTS program_options
    REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries("MultiCompMesher" LINK_PUBLIC ${Boost_LIBRARIES})
    install(TARGETS MultiCompMesher RUNTIME DESTINATION bin)
  else()
    message(STATUS "This program requires the Boost library, and will not be compiled.")
  endif()

  if(ACTIVATE_CONCURRENCY)
    find_package(TBB QUIET)

    if(TBB_FOUND)
      include(CGAL_TBB_support)
      target_link_libraries(MultiCompMesher PRIVATE CGAL::TBB_support)
      add_definitions(-DCGAL_CONCURRENT_MESH_3)
      message(STATUS "NOTICE: Concurrency enabled.")

    else()
      message(STATUS "NOTICE: Unable to find TBB, concurrency will be disabled.")
    endif()
  endif()
else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
