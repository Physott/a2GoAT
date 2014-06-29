cmake_minimum_required (VERSION 2.6)
project(a2GoAT)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()
    
include(cmake/settings.cmake)
include(cmake/doxygen.cmake)

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# require a fairly recent ROOT version
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT "5.30" REQUIRED)

# since all subprojects need ROOT, set that up here
include_directories(${ROOT_INCLUDES})
link_directories(${ROOT_LIBRARY_DIR})

include_directories(inc)
set(ROOT_LIBRARIES "${ROOT_LIBRARIES} -lEG")
message(${ROOT_LIBRARIES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")


add_executable(goat
   inc/GTree.h
   inc/GTreeRawEvent.h
   inc/GTreeDetectorHits.h
   inc/GTreeTagger.h
   inc/GTreeLinPol.h
   inc/GTreeScaler.h
   inc/GTreeParticle.h
   inc/GTreeMeson.h
   inc/GTreeTrigger.h
   inc/GTreeEventParameters.h
   inc/GTreeManager.h
   inc/GConfigFile.h
   inc/GDataChecks.h
   inc/GParticleReconstruction.h
   inc/GMesonReconstruction.h
   inc/GSort.h
   inc/GoAT.h
   src/GTree.cc
   src/GTreeRawEvent.cc
   src/GTreeDetectorHits.cc
   src/GTreeTagger.cc
   src/GTreeLinPol.cc
   src/GTreeScaler.cc
   src/GTreeParticle.cc
   src/GTreeMeson.cc
   src/GTreeTrigger.cc
   src/GTreeEventParameters.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
   src/GDataChecks.cc
   src/GParticleReconstruction.cc
   src/GMesonReconstruction.cc
   src/GSort.cc
   src/GoAT.cc
   src/goat_main.cc
   )
target_link_libraries(goat ${ROOT_LIBRARIES})

add_executable(etap
   inc/GTree.h
   inc/GTreeRawEvent.h
   inc/GTreeDetectorHits.h
   inc/GTreeTagger.h
   inc/GTreeLinPol.h
   inc/GTreeScaler.h
   inc/GTreeParticle.h
   inc/GTreeMeson.h
   inc/GTreeTrigger.h
   inc/GTreeEventParameters.h
   inc/GTreeManager.h
   inc/GConfigFile.h
   inc/GKinFitter.h
   inc/GKinFitterParticle.h
   inc/PHist.h
   inc/PHistEvent.h
   inc/PPhysics.h
   inc/PProtonCheck.h
   inc/P3Meson.h
   inc/MyPhysics.h
   src/GTree.cc
   src/GTreeRawEvent.cc
   src/GTreeDetectorHits.cc
   src/GTreeTagger.cc
   src/GTreeLinPol.cc
   src/GTreeScaler.cc
   src/GTreeParticle.cc
   src/GTreeMeson.cc
   src/GTreeTrigger.cc
   src/GTreeEventParameters.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
   src/GKinFitter.cc
   src/GKinFitterParticle.cc
   src/PHist.cc
   src/PHistEvent.cc
   src/PPhysics.cc
   src/PProtonCheck.cc
   src/P3Meson.cc
   src/MyPhysics.cc
   )
target_link_libraries(etap ${ROOT_LIBRARIES})

