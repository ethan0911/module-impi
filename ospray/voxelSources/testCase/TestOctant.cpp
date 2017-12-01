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

namespace ospray {
  namespace impi { 
    namespace testCase {

      TestOctant::TestOctant(){
      }

      void TestOctant::initOctant(size_t octNum,vec3f* octVertex, float* octWidth,float* octValue)
      {
        std::cout << "Start to Init Octant Value" << std::endl;
        this->octNum = octNum;
        this->octVtxBuffer   = octVertex;
        this->octWidthBuffer = octWidth;
        this->octValueBuffer = octValue;

        // for (size_t i = 0; i < this->octNum; i++) {
        //   Range range;
        //   for (size_t j = 0; j < 8; j++) {
        //     size_t idx = i * 8 + j;
        //     range.extend(this->octValueBuffer[idx]);
        //   }
        //   this->octRange.push_back(range);
        // }

        std::cout << "Done Init Octant Value!" << std::endl;
      }

      TestOctant::~TestOctant(){
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
        float clipping = 250.0f;
        activeVoxels.clear();
        std::cout<<"Filter---------------------------"<<std::endl;
        for (size_t i = 0; i < this->octNum; i++) {
            Range range;
            for (size_t j = 0; j < 8; j++) {
              size_t idx = i * 8 + j;
              range.extend(this->octValueBuffer[idx]);
            }

          auto box = box3fa(this->octVtxBuffer[i],
                              this->octVtxBuffer[i] + vec3f(this->octWidthBuffer[i]));
          //PRINT(box);
          if (range.contains(isoValue) /*&& (box.upper.x < clipping)*/) {
            activeVoxels.push_back(i);
          }
            
        }
        std::cout<<"IsoValue = "<<isoValue<<", ActiveVoxels ="<<activeVoxels.size()<<std::endl;
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        // const Octant& oct = octants[(const int)voxelRef];
        // return box3fa(oct.bounds.lower,oct.bounds.upper);
        return box3fa(this->octVtxBuffer[voxelRef],
                              this->octVtxBuffer[voxelRef] + vec3f(this->octWidthBuffer[voxelRef]));
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        // const Octant& oct = octants[(const int)voxelRef];
        // voxel.bounds = box3fa(oct.bounds.lower,oct.bounds.upper);
        // array3D::for_each(vec3i(2),[&](const vec3i vtx){
        //   voxel.vtx[vtx.z][vtx.y][vtx.x] =
        //   oct.vertexValue[vtx.z][vtx.y][vtx.x];
        // });
        //size_t octID = (size_t)VoxelRef;
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
