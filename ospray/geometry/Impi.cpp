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

#include "Impi.h"
// 'export'ed functions from the ispc file:
#include "Impi_ispc.h"
// ospray core:
#include <ospray/common/Data.h>

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

    /*! constructor - will create the 'ispc equivalent' */
    Impi::Impi()
    {
      /*! create the 'ispc equivalent': ie, the ispc-side class that
        implements all the ispc-side code for intersection,
        postintersect, etc. See Impi.ispc */
      this->ispcEquivalent = ispc::Impi_create(this);

      // note we do _not_ yet do anything else here - the actual input
      // data isn't available to use until 'commit()' gets called
    }

    /*! destructor - supposed to clean up all alloced memory */
    Impi::~Impi()
    {
      ispc::Impi_destroy(ispcEquivalent);
    }

    /*! commit - this is the function that parses all the parameters
      that the app has proivded for this geometry. In this simple
      example we're looking for a single parameter named 'patches',
      which is supposed to contain a data array of all the patches'
      control points */
    void Impi::commit()
    {
      // this->voxelData = getParamData("voxel");

      /* assert that some valid input data is available */
    }


    extern "C" void externC_getCell(Cell &cell,
                                    const LogicalVolume *volume,
                                    const vec3i &cellIdx)
    {
      volume->getCell(cell,cellIdx);
    }

    /*! 'finalize' is what ospray calls when everything is set and
        done, and a actual user geometry has to be built */
    void Impi::finalize(Model *model)
    {
      // sanity check if a patches data was actually set!
      float isoValue = 20.f;//0.4f;

#if 1
      std::cout << "loading test data-set, and testing generation of iso-voxels" << std::endl;
      volume = VolumeT<float>::loadRAW("density_064_064_2.0.raw",vec3i(64));
      seg    = VolumeT<float>::loadRAW("density_064_064_2.0_seg.raw",vec3i(64));

      volume->filterVoxels(hotCells,[&](const LogicalVolume *v, const vec3i &idx) {
          return
            (volume->getCell(idx).getRange().contains(isoValue)
             &&
             seg->getCell(idx).getRange().contains(128)
             );
        });
      
      std::cout << "asking ISPC to build a bvh over the hot cells..." << std::endl;
      vec3i dims = volume->getDims();
      ispc::Impi_finalize_embreeBVHoverHotCells(getIE(),model->getIE(),
                                                (uint64_t*)&hotCells[0],hotCells.size(),
                                                (ispc::vec3i &)dims,
                                                volume.get(),
                                                isoValue);
#else
      /* get the acual 'raw' pointer to the data (ispc doesn't konw
         what to do with the 'Data' abstraction calss */
      void *voxelDataPointer = voxelData->data;
      ispc::Impi_finalize_testCell(getIE(),model->getIE(),
                                   (float*)voxelDataPointer,
                                   isoValue);
#endif
    }


    /*! maybe one of the most important parts of this example: this
        macro 'registers' the Impi class under the ospray
        geometry type name of 'bilinear_patches'.

        It is _this_ name that one can now (assuming the module has
        been loaded with ospLoadModule(), of course) create geometries
        with; i.e.,

        OSPGeometry geom = ospNewGeometry("bilinear_patches") ;
    */
    OSP_REGISTER_GEOMETRY(Impi,impi);

  } // ::ospray::bilinearPatch
} // ::ospray
