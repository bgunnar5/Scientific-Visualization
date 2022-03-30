cmake_minimum_required(VERSION 3.10)

PROJECT(isosurfaces)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)

add_executable(isosurfaces MACOSX_BUNDLE isosurfaces)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(isosurfaces ${VTK_LIBRARIES})
else()
  target_link_libraries(isosurfaces vtkHybrid)
endif()


