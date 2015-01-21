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
   inc/GHistManager.h
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
   src/GHistManager.cc
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


add_executable(goat_detailed
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
   inc/GHistManager.h
   inc/GTreeManager.h
   inc/GConfigFile.h
   inc/GDataChecks.h
   inc/GParticleReconstruction.h
   inc/GMesonReconstruction_6and7gamma.h
   inc/GSort.h
   inc/GoAT_Detailed.h
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
   src/GHistManager.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
   src/GDataChecks.cc
   src/GParticleReconstruction.cc
   src/GMesonReconstruction_6and7gamma.cc
   src/GSort.cc
   src/GoAT_Detailed.cc
   src/goat_detailed_main.cc
   )
target_link_libraries(goat_detailed ${ROOT_LIBRARIES})


add_executable(pi0-example
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
   inc/GHistManager.h
   inc/GTreeManager.h
   inc/GConfigFile.h
   inc/GHistScaCor.h
   inc/GHistScaCor2.h
   inc/GHistScaCor3.h
   inc/GHistBGSub.h
   inc/GHistBGSub2.h
   inc/GHistBGSub3.h
   inc/GHistTaggerBinning.h
   inc/GHistTaggerBinning2.h
   inc/GH1.h
   inc/PPhysics.h
   inc/PPi0Example.h
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
   src/GHistManager.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
   src/GHistScaCor.cc
   src/GHistScaCor2.cc
   src/GHistScaCor3.cc
   src/GHistBGSub.cc
   src/GHistBGSub2.cc
   src/GHistBGSub3.cc
   src/GHistTaggerBinning.cc
   src/GHistTaggerBinning2.cc
   src/PPhysics.cc
   src/PPi0Example.cc
   src/pi0_example_main.cc
   )
target_link_libraries(pi0-example ${ROOT_LIBRARIES})


add_executable(MyPhysics
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
   inc/GHistManager.h
   inc/GTreeManager.h
   inc/GConfigFile.h
   inc/GHistScaCor.h
   inc/GHistScaCor2.h
   inc/GHistScaCor3.h
   inc/GHistBGSub.h
   inc/GHistBGSub2.h
   inc/GHistBGSub3.h
   inc/GHistTaggerBinning.h
   inc/GHistTaggerBinning2.h
   inc/GH1.h
   inc/GHistFit.h
   inc/GKinFitterConverter.h
   inc/GKinFitter.h
   inc/GHistEvent.h
   inc/GCheckProton.h
   inc/GAnalysis3Mesons.h
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
   src/GHistManager.cc
   src/GTreeManager.cc
   src/GConfigFile.cc
   src/GHistScaCor.cc
   src/GHistScaCor2.cc
   src/GHistScaCor3.cc
   src/GHistBGSub.cc
   src/GHistBGSub2.cc
   src/GHistBGSub3.cc
   src/GHistTaggerBinning.cc
   src/GHistTaggerBinning2.cc
   src/GHistFit.cc
   src/GKinFitterConverter.cc
   src/GKinFitter.cc
   src/GHistEvent.cc
   src/GCheckProton.cc
   src/GAnalysis3Mesons.cc
   src/MyPhysics.cc
   src/MyPhysics_main.cc
   )
target_link_libraries(MyPhysics ${ROOT_LIBRARIES})

# install some scripts to the bin directory
# by creating symlinks
file(GLOB CORE_EXTRA_SCRIPTS "scripts/*")
foreach(f ${CORE_EXTRA_SCRIPTS})
  get_filename_component(f_name ${f} NAME)
  add_custom_target(link_${f_name} ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${f}" "${EXECUTABLE_OUTPUT_PATH}/${f_name}")
endforeach()
