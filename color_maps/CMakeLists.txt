cmake_minimum_required(VERSION 3.10)

PROJECT(color_maps)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)

add_executable(color_maps MACOSX_BUNDLE color_maps)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(color_maps ${VTK_LIBRARIES})
else()
  target_link_libraries(color_maps vtkHybrid)
endif()


