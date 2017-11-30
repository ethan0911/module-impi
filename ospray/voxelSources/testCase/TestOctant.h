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
#include "ospcommon/range.h"
#include "ospcommon/array3D/for_each.h"
#include "../..//geometry/Impi.h"

namespace ospray {
  namespace impi { 
    namespace testCase {

      typedef ospcommon::range_t<float> Range;
      
      struct Octant
      {
        box3f bounds;
        float width;
        float vertexValue[2][2][2];

        inline Range getRange() const {
          Range range;
          array3D::for_each(vec3i(2),[&](const vec3i idx) {
              range.extend(vertexValue[idx.z][idx.y][idx.x]);
            });
          return range;
        }
      };

      /*! implements a simple (vertex-cenetred) AMR test case
          consisting of a 2x2x2-cell base level in which one of the
          cells is refined into another 2x2x2-cell finer level */
      struct TestOctant: public Impi::VoxelSource {

        TestOctant();
        ~TestOctant();

        //void parseOctant(std::string fileName);
        /*! create lits of *all* voxel (refs) we want to be considered for interesction */
        virtual void   getActiveVoxels(std::vector<VoxelRef> &activeVoxels, float isoValue) const override;
        
        /*! compute world-space bounds for given voxel */
        virtual box3fa getVoxelBounds(const VoxelRef voxelRef) const override;
        
        /*! get full voxel - bounds and vertex values - for given voxel */
        virtual Impi::Voxel  getVoxel(const VoxelRef voxelRef) const override;

        void initOctant(size_t octNum,vec3f* octVertex, float* octValue);

        std::vector<Octant> octants;

        size_t octNum;
        vec3f* octVtxBuffer;
        float* octValueBuffer;
        std::vector<Range> octRange;
      };

    }
  }
}
