// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

/*! \file ospray/geometry/Impi.h Defines a new ospray
  geometry type (of name 'bilinear_patches'). Input to the geometry is
  a single data array (named 'patches') that consists of for vec3fs
  per patch. */

// ospcomon: vec3f, box3f, etcpp - generic helper stuff
#include <ospcommon/vec.h>
#include <ospcommon/box.h>
// ospray: everything that's related to the ospray ray tracing core
#include <ospray/geometry/Geometry.h>
#include <ospray/common/Model.h>

// OUR includes
#include "../common/Volume.h"

/*! _everything_ in the ospray core universe should _always_ be in the
  'ospray' namespace. */
namespace ospray {

  /*! though not required, it is good practice to put any module into
    its own namespace (isnide of ospray:: ). Unlike for the naming of
    library and init function, the naming for this namespace doesn't
    particularlly matter. E.g., 'bilinearPatch', 'module_blp',
    'bilinar_patch' etc would all work equally well. */
  namespace impi {
    // import ospcommon component - vec3f etc
    using namespace ospcommon;

    /*! a geometry type that implements (a set of) bi-linear
      patches. This implements a new ospray geometry, and as such has
      to

      a) derive from ospray::Geometry
      b) implement a 'commit()' message that parses the
         parameters/data arrays that the app has specified as inputs
      c) create an actual ospray geometry instance with the
         proper intersect() and postIntersect() functions.

      Note that how this class is called does not particularly matter;
      all that matters is under which name it is registered in the cpp
      file (see comments on OSPRAY_REGISTER_GEOMETRY)
    */
    struct Impi : public ospray::Geometry
    {
      /*! constructor - will create the 'ispc equivalent' */
      Impi();

      /*! destructor - supposed to clean up all alloced memory */
      virtual ~Impi() override;

      /*! the commit() message that gets called upon the app calling
          "ospCommit(<thisGeometry>)" */
      virtual void commit() override;

      /*! 'finalize' is what ospray calls when everything is set and
        done, and a actual user geometry has to be built */
      virtual void finalize(Model *model) override;

      /*! for the case where we build an embree bvh over all hot
          cells, this is the vector that stores them.. */
      std::vector<CellRef> hotCells;
      
      std::shared_ptr<LogicalVolume> volume;
    };

  } // ::ospray::bilinearPatch
} // ::ospray

