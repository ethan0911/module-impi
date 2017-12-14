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

#include "TestOctant.h"
#include "ospcommon/tasking/parallel_for.h"
#include "time.h"



#ifndef speedtest__
#define speedtest__(data)                                           \
  for (long blockTime = NULL;                                       \
       (blockTime == NULL ? (blockTime = clock()) != NULL : false); \
       printf("Calculate Time: %.9fs \n", (double)(clock() - blockTime) / CLOCKS_PER_SEC))
#endif

namespace ospray {
  namespace impi { 
    namespace testCase {

      TestOctant::TestOctant(){
      }

      void TestOctant::initOctant(ospray::AMRVolume * amrDataNode)
      {
        assert(amrDataNode != NULL);
        std::cout << "Start to Init Octant Value" << std::endl;
        this->octNum         = amrDataNode->accel->octNum;
        this->octVtxBuffer   = (vec3f*)amrDataNode->accel->octVertices.data();
        this->octWidthBuffer = (float *)amrDataNode->accel->octWidth.data();
        this->octValueBuffer = amrDataNode->accel->octVerticeValue;

        // this->clappingBox = box3fa(amrDataNode->accel->worldBounds.lower,
        //                           0.5f * amrDataNode->accel->worldBounds.upper);

        auto wb = amrDataNode->accel->worldBounds;
        PRINT(wb);
        vec3f center = wb.center();
        box3fa b1 = box3fa(wb.lower,wb.upper);
        b1.upper.x *= 0.5f;
        box3fa b2 = box3fa(vec3f(center.x,0.f,0.f),vec3f(wb.upper.x,center.y,center.z));
        box3fa b3 = box3fa(vec3f(center.x,center.y,0.f),wb.upper);
        this->clapBoxes.push_back(b1);
        this->clapBoxes.push_back(b2);
        this->clapBoxes.push_back(b3);

        speedtest__("Speed: ")
        {
          this->octRange.resize(this->octNum);
          for (size_t i = 0; i < this->octNum; i++) {
            getOctrange(i, &this->octRange[i]);
          }
        }
        std::cout << "Done Init Octant Value!" << std::endl;
      }

      TestOctant::~TestOctant()
      {
        if (octVtxBuffer != NULL)
          delete[] octVtxBuffer;
        if(octWidthBuffer != NULL)
          delete[] octWidthBuffer;
        if (octValueBuffer != NULL)
          delete[] octValueBuffer;
      }

      /*! create lits of *all* voxel (refs) we want to be considered for
       * interesction */
      void TestOctant::getActiveVoxels(std::vector<VoxelRef> &activeVoxels, float isoValue) const 
      {
        activeVoxels.clear();
        std::cout<<"Filter---------------------------"<<std::endl;

        for (size_t i = 0; i < this->octNum; i++) {
          auto box = box3fa(this->octVtxBuffer[i],this->octVtxBuffer[i] + vec3f(this->octWidthBuffer[i]));
          if (octRange[i].contains(isoValue) /*&& isInClapBox(box)*/) {
            activeVoxels.push_back(i);
          }
        }
   
        std::cout<<"IsoValue = "<<isoValue<<", ActiveVoxels ="<<activeVoxels.size()<<std::endl;
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        return box3fa(this->octVtxBuffer[voxelRef],
                              this->octVtxBuffer[voxelRef] + vec3f(this->octWidthBuffer[voxelRef]));
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        size_t startIdx = voxelRef * 8;
        voxel.bounds = box3fa(this->octVtxBuffer[voxelRef],
                              this->octVtxBuffer[voxelRef] + vec3f(this->octWidthBuffer[voxelRef]));
        voxel.vtx[0][0][0] = this->octValueBuffer[startIdx++];
        voxel.vtx[0][0][1] = this->octValueBuffer[startIdx++];
        voxel.vtx[0][1][0] = this->octValueBuffer[startIdx++];
        voxel.vtx[0][1][1] = this->octValueBuffer[startIdx++];
        voxel.vtx[1][0][0] = this->octValueBuffer[startIdx++];
        voxel.vtx[1][0][1] = this->octValueBuffer[startIdx++];
        voxel.vtx[1][1][0] = this->octValueBuffer[startIdx++];
        voxel.vtx[1][1][1] = this->octValueBuffer[startIdx++];
        return voxel;
      }
    }
  }
}
