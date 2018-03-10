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
if (OSPRAY_MODULE_IMPLICIT_ISOSURFACES_VIEWER)
  ## ========================================================================== ##
  ## 
  ## ========================================================================== ##
  find_package(OpenGL REQUIRED)
  ospray_create_application(ospImplicitIsoSurfaceViewer
    viewer/impiViewer.cpp
    LINK
    ospray
    ospray_common
    ospray_imgui3d_sg
    ${GLUT_LIBRARIES}
    ${OPENGL_LIBRARIES})
endif(OSPRAY_MODULE_IMPLICIT_ISOSURFACES_VIEWER)

if (OSPRAY_MODULE_IMPLICIT_ISOSURFACES_BENCHMARKER)

  set(USE_VIEWER 0)
  if (OSPRAY_MODULE_IMPLICIT_ISOSURFACES_VIEWER)
    if (TFN_MODULE_ROOT)
      ## ====================================================================== ##
      ##
      ## ====================================================================== ##
      find_package(OpenGL REQUIRED)
      find_package(GLFW REQUIRED)
      set(USE_VIEWER 1)
      set(VIEWER_INC
	bench/opengl/external
	bench/opengl/external/ImGui/
	bench/opengl/external/ImGui/examples/opengl3_example/
	${TFN_MODULE_ROOT}/src)
      set(VIEWER_SRC
	bench/opengl/viewer.cpp
	bench/opengl/external/glad.c
	bench/opengl/external/ImGui/imgui.cpp
	bench/opengl/external/ImGui/imgui_draw.cpp
	bench/opengl/external/ImGui/imgui_demo.cpp
	bench/opengl/external/ImGui/examples/opengl3_example/imgui_impl_glfw_gl3.cpp
	${TFN_MODULE_ROOT}/src/TransferFunctionModule.cpp
	${TFN_MODULE_ROOT}/src/common/jsoncpp.cpp
	${TFN_MODULE_ROOT}/src/widgets/TransferFunctionWidget.cpp)
      set(VIEWER_LIB 
	${GLFW_LIBRARIES}
	${GLUT_LIBRARIES}
	${OPENGL_LIBRARIES})
    endif (TFN_MODULE_ROOT)
  endif (OSPRAY_MODULE_IMPLICIT_ISOSURFACES_VIEWER)

  ## ======================================================================== ##
  ## 
  ## ======================================================================== ##
  include_directories(${VIEWER_INC})
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
  endif()
  ospray_create_application(ospImplicitIsoSurfaceBench
    bench/impiBench.cpp
    bench/impiReader.cpp
    bench/loader/meshloader.cpp
    ${VIEWER_SRC}
    LINK
    ospray
    ospray_common
    ospray_xml
    ${VIEWER_LIB}
    ${HDF5_LIBRARIES})
  set_target_properties(ospImplicitIsoSurfaceBench
    PROPERTIES 
    CXX_STANDARD 14
    COMPILE_DEFINITIONS
    USE_VIEWER=${USE_VIEWER})

endif (OSPRAY_MODULE_IMPLICIT_ISOSURFACES_BENCHMARKER)
