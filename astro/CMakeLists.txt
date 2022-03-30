cmake_minimum_required(VERSION 3.10)

PROJECT(astro)
SET(VTK_DIR ADD PATH TO VTK BUILD HERE)

find_package(VTK REQUIRED)

add_executable(astro MACOSX_BUNDLE astro)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(astro ${VTK_LIBRARIES})
else()
  target_link_libraries(astro vtkHybrid)
endif()


