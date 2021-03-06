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

# Define an option to enable/disable this module.
#
# Though not actually required we guard each module 'moduleName' with
# a OSPRAY_MODULE_MODULENAME flag that enables/disables its building.
#
OPTION(OSPRAY_MODULE_IMPLICIT_ISOSURFACES
  "Build module the provides implicit iso-surfaces in structured grids" ON)
OPTION(OSPRAY_MODULE_IMPI_VIEWER 
  "Build the ImGui/OpenGL viewer implicit iso-surfaces" OFF)
OPTION(OSPRAY_MODULE_IMPI_BENCH_WIDGET
  "Build a simple OpenGL viewer for implicit iso-surfaces" ON)
OPTION(OSPRAY_MODULE_IMPI_BENCH_MARKER
  "Build the benchmarker for implicit iso-surfaces" ON)

INCLUDE_DIRECTORIES_ISPC(${EMBREE_INCLUDE_DIRS}/embree3)
INCLUDE_DIRECTORIES(${EMBREE_INCLUDE_DIRS}/embree3)


IF (OSPRAY_MODULE_IMPLICIT_ISOSURFACES)
  ADD_SUBDIRECTORY(ospray) 
  ADD_SUBDIRECTORY(apps)
ENDIF ()
