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
include_directories(${ROOT_INCLUDES} ${CMAKE_SOURCE_DIR}/inc)
link_directories(${ROOT_LIBRARY_DIR})

include_directories(inc)
set(ROOT_LIBRARIES "${ROOT_LIBRARIES} -lEG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall")

set(GOAT_BASE
   inc/GTree.h
   inc/GTreeTrack.h
   inc/GTreeDetectorHits.h
   inc/GTreeTagger.h
   inc/GTreeLinPol.h
   inc/GTreeScaler.h
   inc/GTreeParticle.h
   inc/GTreeMeson.h
   inc/GTreeTrigger.h
   inc/GTreeSetupParameters.h
   inc/GTreeEventParameters.h
   inc/GTreeA2Geant.h
   inc/GHistManager.h
   inc/GTreeManager.h
   inc/GConfigFile.h   
   src/GTree.cc
   src/GTreeTrack.cc
   src/GTreeDetectorHits.cc
   src/GTreeTagger.cc
   src/GTreeLinPol.cc
   src/GTreeScaler.cc
   src/GTreeParticle.cc
   src/GTreeMeson.cc
   src/GTreeTrigger.cc
   src/GTreeSetupParameters.cc
   src/GTreeEventParameters.cc
   src/GTreeA2Geant.cc
   src/GHistManager.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
)

set(GOAT_GHIST
   inc/GH1.h
   inc/GHistScaCor.h
   src/GHistScaCor.cc
   inc/GHistScaCor2.h
   src/GHistScaCor2.cc
   inc/GHistScaCor3.h
   src/GHistScaCor3.cc
   inc/GHistBGSub.h
   src/GHistBGSub.cc
   inc/GHistBGSub2.h
   src/GHistBGSub2.cc
   inc/GHistBGSub3.h
   src/GHistBGSub3.cc
   inc/GHistTaggerBinning.h
   src/GHistTaggerBinning.cc
   inc/GHistTaggerBinning2.h
   src/GHistTaggerBinning2.cc
)

if( PLUTO )

    find_package(Pluto)

    if( PLUTO_FOUND )

           set(GOAT_BASE
              ${GOAT_BASE}
              inc/GTreePluto.h
              src/GTreePluto.cc
           )

           add_definitions( -DhasPluto )

           include_directories(${PLUTO_INCLUDE_PATH})
           link_directories(${PLUTO_LIBRARY_PATH})

           set(LIBS ${LIB} ${PLUTO_LIBRARY})
    endif()

else()
    message(STATUS "Pluto disabled")
endif()


add_executable(goat
   ${GOAT_BASE}
   inc/GDataChecks.h
   src/GDataChecks.cc
   inc/GParticleReconstruction.h
   src/GParticleReconstruction.cc
   inc/GMesonReconstruction.h
   src/GMesonReconstruction.cc
   inc/GSort.h
   src/GSort.cc
   inc/GoAT.h
   src/GoAT.cc
   src/goat_main.cc
   )
target_link_libraries(goat ${LIBS} ${ROOT_LIBRARIES})

add_executable(pi0-example
   ${GOAT_BASE}
   ${GOAT_GHIST}
   inc/PPhysics.h
   src/PPhysics.cc
   inc/PPi0Example.h
   src/PPi0Example.cc
   src/pi0_example_main.cc
   )
target_link_libraries(pi0-example ${LIBS} ${ROOT_LIBRARIES})

add_executable(gh-example
   ${GOAT_BASE}
   ${GOAT_GHIST}
   inc/GHExample.h
   src/GHExample.cc
   src/gh_example_main.cc
   )
target_link_libraries(gh-example ${LIBS} ${ROOT_LIBRARIES})

# install some scripts to the bin directory
# by creating symlinks
file(GLOB CORE_EXTRA_SCRIPTS "scripts/*")
foreach(f ${CORE_EXTRA_SCRIPTS})
  get_filename_component(f_name ${f} NAME)
  add_custom_target(link_${f_name} ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${f}" "${EXECUTABLE_OUTPUT_PATH}/${f_name}")
endforeach()
