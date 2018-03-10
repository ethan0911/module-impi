## ======================================================================== ##
## Copyright 2009-2017 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

ospray_disable_compiler_warnings()

# this subdirectory demonstrates on how to _use_ demo module from
# within a simple viewer

# the actual example viewer. Note this viewer does _NOT_ link to the
# actual module library - the module should (and _does_, in this
# example) get loaded through ospLoadModule(...), which will make sure
# the module gets properly oaded on all worker nodes, even if those
# live on different nodes, accelerator cards, memory spaces, etc.
## ======================================================================== ##
## SG viewer
## ======================================================================== ##
if (OSPRAY_MODULE_IMPI_VIEWER)
  find_package(OpenGL REQUIRED)
  ospray_create_application(ospImplicitIsoSurfaceViewer
    viewer/impiViewer.cpp
    LINK
    ospray
    ospray_common
    ospray_imgui3d_sg
    ${OPENGL_LIBRARIES})
endif(OSPRAY_MODULE_IMPI_VIEWER)

## ====================================================================== ##
##
## ====================================================================== ##
if (OSPRAY_MODULE_IMPI_BENCH_WIDGET)
  find_package(OpenGL REQUIRED)
  find_package(GLFW REQUIRED)
  include_directories(
    bench/opengl/external
    bench/opengl/external/ImGui
    bench/opengl/external/ImGui/examples/opengl3_example
    ${PROJECT_SOURCE_DIR}/apps/exampleViewer/widgets/TransferFunctionModule/src
    ${GLFW_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR})
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
  endif()
  ospray_create_application(ospImplicitIsoSurfaceWidget
    bench/impiBench.cpp
    bench/impiReader.cpp
    bench/loader/meshloader.cpp
    bench/opengl/viewer.cpp
    bench/opengl/external/glad.c
    bench/opengl/external/ImGui/imgui.cpp
    bench/opengl/external/ImGui/imgui_draw.cpp
    bench/opengl/external/ImGui/imgui_demo.cpp
    bench/opengl/external/ImGui/examples/opengl3_example/imgui_impl_glfw_gl3.cpp
    ${PROJECT_SOURCE_DIR}/apps/exampleViewer/widgets/TransferFunctionModule/src/TransferFunctionModule.cpp
    ${PROJECT_SOURCE_DIR}/apps/exampleViewer/widgets/TransferFunctionModule/src/common/jsoncpp.cpp
    ${PROJECT_SOURCE_DIR}/apps/exampleViewer/widgets/TransferFunctionModule/src/widgets/TransferFunctionWidget.cpp
    LINK
    ospray
    ospray_common
    ospray_xml
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${HDF5_LIBRARIES})
  set_target_properties(ospImplicitIsoSurfaceWidget
    PROPERTIES 
    CXX_STANDARD 14
    COMPILE_DEFINITIONS
    USE_VIEWER=1)
endif (OSPRAY_MODULE_IMPI_BENCH_WIDGET)

## ==================================================================== ##
## Benchmarker
## ==================================================================== ##
if (OSPRAY_MODULE_IMPI_BENCH_MARKER)
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
  endif()
  ospray_create_application(ospImplicitIsoSurfaceBench
    bench/impiBench.cpp
    bench/impiReader.cpp
    bench/loader/meshloader.cpp
    LINK
    ospray
    ospray_common
    ospray_xml
    ${HDF5_LIBRARIES})
  set_target_properties(ospImplicitIsoSurfaceBench
    PROPERTIES 
    CXX_STANDARD 14
    COMPILE_DEFINITIONS
    USE_VIEWER=0)
endif (OSPRAY_MODULE_IMPI_BENCH_MARKER)
