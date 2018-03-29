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
#include "volume/amr/AMRVolume.h"
#include "volume/amr/AMRAccel.h"
#include "../../geometry/Impi.h"
#include <limits>

namespace ospray {
  namespace impi { 
    namespace testCase {

      typedef ospcommon::range_t<float> Range;
      
      struct Octant
      {	
	vec3f lowerleft;
	float cellwidth;
        float vtx[2][2][2];
      };

      /*! implements a simple (vertex-cenetred) AMR test case
          consisting of a 2x2x2-cell base level in which one of the
          cells is refined infilterActiveVoxelsto another 2x2x2-cell
	  finer level */
      struct TestOctant : public Impi::VoxelSource
      {
      public:	
	/*! constructors and distroctors */
        TestOctant();
        virtual ~TestOctant();

        /*! compute world-space bounds for given voxel */
        virtual box3fa getVoxelBounds(const VoxelRef voxelRef) const override;

        /*! get full voxel - bounds and vertex values - for given voxel */
        virtual Impi::Voxel getVoxel(const VoxelRef voxelRef) const override;
        	
	
	/*! compute active voxels (called in Impi.cpp file) */
        virtual void getActiveVoxels(std::vector<VoxelRef> &activeVoxels, 
				     float isoValue) const override;

      public: 

        void buildActiveVoxels(std::vector<VoxelRef> &activeVoxels, 
			       float isoValue);


	/*! initialization */
        void initOctant(ospray::AMRVolume *amr);
	


	/*! check if the voxel is inside the clip box */
	bool inClipBox(const box3f& box) const {
	  return inClipBox(box3fa(box.lower, box.upper));
	}
	bool inClipBox(const box3fa& box) const {
          bool result = false;
          for(auto clipbox: clipBoxes) {
            if (touchingOrOverlapping(clipbox, box)) {
              result = true;
              break;
            }
          }
          return result;
        }

      private:
        std::vector<Octant> octants; // active octant list
        std::vector<box3fa> clipBoxes;
	const AMRVolume *amrVolumePtr = nullptr;


      private:
	//TODO things below should be cleaned
	
        //void parseOctant(std::string fileName);
        /*! create lits of *all* voxel (refs) we want to be considered for 
	  interesction */
        size_t octNum = 0;

        void initOctantValue(ospray::AMRVolume *amr);



        void buildOctant(ospray::AMRVolume *amr);
        void buildOctantByLeaf(std::unique_ptr<ospray::amr::AMRAccel> &accel,
                               ospray::amr::AMRAccel::Leaf &lf,
                               std::vector<vec3f> *outOctLowV,
                               std::vector<float> *outOctW);



        std::vector<vec3f> octVertices;
        std::vector<float> octWidth;
        float *octValueBuffer;

        std::vector<Range> octRange;


	
	void getOctrange(size_t octID,Range* range)
        {
          for (size_t i = 0; i < 8; i++) {
            size_t idx = octID * 8 + i;
            range->extend(this->octValueBuffer[idx]);
          }
        }
      };

    }  // namespace testCase
  }
}
