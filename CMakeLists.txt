cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name CaTS)
project(${name})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
add_definitions(-DG4ANALYSIS_USE) 
#find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "CaTS: ROOT package not found. --> CaTS disabled")
  return()
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

REFLEX_GENERATE_DICTIONARY(Classes include/Classes.hh
                           SELECTION xml/selection.xml
                           OPTIONS --deep --iocomments --fail_on_warnings)

add_library(ClassesDict SHARED ${sources} Classes.cpp)

set(libsuffix .so)

set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES}
      SUFFIX ${libsuffix})

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set_target_properties(ClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

SET(ROOT_LIBRARIES "${ROOT_LIBRARIES} -lGui")

target_link_libraries(ClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex Reflex)

add_executable(${name} ${name}.cc )
target_link_libraries(${name} ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

add_executable(readHits readHits.cc )
target_link_libraries(readHits ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(drawShower drawShower.cc )
#target_link_libraries(drawShower ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

add_executable(drawShowerGUI drawShowerGUI.cc )
target_link_libraries(drawShowerGUI ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(DumpHeader DumpHeader.cc )
#target_link_libraries(DumpHeader ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(uniformity uniformity.cc )
#target_link_libraries(uniformity ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(profile profile.cc )
#target_link_libraries(profile ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)
 
#add_executable(Calibration Calibration.cc )
#target_link_libraries(Calibration ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(Hits Hits.cc )
#target_link_libraries(Hits ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(ialength ialength.cc )
#target_link_libraries(ialength ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

#add_executable(response response.cc )
#target_link_libraries(response ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

install (TARGETS ${name} readHits drawShowerGUI DESTINATION bin)
install (TARGETS ${name} ClassesDict DESTINATION lib)
install (FILES vis.mac pythia_event.data SimpleExtension.xsd DESTINATION bin)
