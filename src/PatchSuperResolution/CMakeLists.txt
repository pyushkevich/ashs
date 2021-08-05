project(SuperResolution)
cmake_minimum_required(VERSION 2.8)


#SET(ITK_DIR "/Users/long/Development/bin/ITK_install/lib/cmake/ITK-4.10" CACHE PATH "ITK directory override" FORCE)
#SET(VTK_DIR "/Users/long/Development/bin/VTK_install/lib/cmake/vtk-7.1" CACHE PATH "VTK directory override" FORCE)
SET(VTK_DIR "/home/cwang/Libraries/lib/cmake/vtk-7.1" CACHE PATH "VTK directory override" FORCE)


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
#if (ITKVtkGlue_LOADED)
#  find_package(VTK REQUIRED)
#  include(${VTK_USE_FILE})
#else()
#  find_package(ItkVtkGlue REQUIRED)
#  include(${ItkVtkGlue_USE_FILE})
#  set(Glue ItkVtkGlue)
#endif()

include_directories (src)
include_directories(./include/)
#include_directories(./include/util)

#SET( SOURCES 
#	${PROJECT_SOURCE_DIR}/src/
#)

ADD_EXECUTABLE(NLMDenoise src/MABONLM3D.cxx)
target_link_libraries(NLMDenoise ${ITK_LIBRARIES})

#target_link_libraries(NLMDenoise
#    ${VTK_LIBRARIES} ${ITK_LIBRARIES})
#target_link_libraries(SuperResolution ${ITK_LIBRARIES})

ADD_EXECUTABLE(NLMUpsample src/NLMUpsample.cxx)
target_link_libraries(NLMUpsample ${ITK_LIBRARIES})
#target_link_libraries(NLMUpsample
#    ${VTK_LIBRARIES} ${ITK_LIBRARIES})
