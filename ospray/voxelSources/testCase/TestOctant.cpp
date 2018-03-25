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
#include "ospcommon/utility/getEnvVar.h"
#include "time.h"
#include <ospray/common/Data.h>
#include "ospray/AMRVolume_ispc.h"
#include "ospray/method_finest_ispc.h"
#include "ospray/method_current_ispc.h"
#include "ospray/method_octant_ispc.h"

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

      void TestOctant::initOctant(ospray::AMRVolume * amr)
      {
        assert(amr != NULL);
        std::cout << "Start to Init Octant Value" << std::endl;
        
        buildOctant(amr);
        initOctantValue(amr);

        auto wb = amr->accel->worldBounds;

        // vec3f center = wb.center();
        // box3fa b1 = box3fa(wb.lower,wb.upper);
        // b1.upper.x *= 0.5f;
        // box3fa b2 = box3fa(vec3f(center.x,0.f,0.f),vec3f(wb.upper.x,center.y,center.z));
        // box3fa b3 = box3fa(vec3f(center.x,center.y,0.f),wb.upper);
        // this->clapBoxes.push_back(b1);
        // this->clapBoxes.push_back(b2);
        // this->clapBoxes.push_back(b3);

        box3fa b1 = box3fa(wb.lower, wb.upper);
        b1.upper.z *= 0.5f;
        this->clapBoxes.push_back(b1);

        speedtest__("Speed: ")
        {
          this->octRange.resize(this->octNum);
          for (size_t i = 0; i < this->octNum; i++) {
            getOctrange(i, &this->octRange[i]);
          }
        }
        std::cout << "Done Init Octant Value!" << std::endl;
      }

      void TestOctant::buildOctantByLeaf(std::unique_ptr<ospray::amr::AMRAccel> &accel,
                                         ospray::amr::AMRAccel::Leaf &lf,
                                         std::vector<vec3f> *outOctLowV,
                                         std::vector<float> *outOctW)
      {
        std::vector<vec3f> oVertice;
        std::vector<float> oWidth;
        float cw = lf.brickList[0]->cellWidth;
        vec3f &lower = lf.bounds.lower;
        vec3f &upper = lf.bounds.upper;
#if 1
        // Octant method vertex array
        const float halfCW = 0.5f * cw;
        float fz, fy, fx;
        // add inner cells
        for (fz = lower.z + halfCW; fz < upper.z - halfCW; fz += cw)
          for (fy = lower.y + halfCW; fy < upper.y - halfCW; fy += cw)
            for (fx = lower.x + halfCW; fx < upper.x - halfCW; fx += cw) {
              oWidth.push_back(cw);
              oVertice.push_back(vec3f(fx, fy, fz));
            }
        //bottom top boundray cells
        for (fy = lower.y; fy < upper.y; fy += halfCW)
          for (fx = lower.x; fx < upper.x; fx += halfCW) {
            fz = lower.z;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
            fz = upper.z - halfCW;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
          }
        // left right boundray cells
        for (fz = lower.z; fz < upper.z; fz += halfCW)
          for (fy = lower.y; fy < upper.y; fy += halfCW) {
            fx = lower.x;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
            fx = upper.x - halfCW;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
          }
        // front back boundary cells
        for (fz = lower.z; fz < upper.z; fz += halfCW)
          for (fx = lower.x; fx < upper.x; fx += halfCW) {
            fy = lower.y;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
            fy = upper.y - halfCW;
            oWidth.push_back(halfCW);
            oVertice.push_back(vec3f(fx, fy, fz));
          }
#elif 0
        //finest VA
        cw = (cw == 1.0f || cw == 0.125f)
                 ? cw
                 : accel->level[accel->level.size()-2].cellWidth;

        for (float fz = lower.z; fz < upper.z; fz += cw)
          for (float fy = lower.y; fy < upper.y; fy += cw)
            for (float fx = lower.x; fx < upper.x; fx += cw) {
              oWidth.push_back(cw);
              oVertice.push_back(vec3f(fx, fy, fz));
            }
#else
      for (float fz = lf.bounds.lower.z; fz < lf.bounds.upper.z;
           fz += 0.5f * cw)
        for (float fy = lf.bounds.lower.y; fy < lf.bounds.upper.y;
             fy += 0.5f * cw)
          for (float fx = lf.bounds.lower.x;
               fx < lf.bounds.upper.x;
               fx += 0.5f * cw) {
            oWidth.push_back(0.5f * cw);
            oVertice.push_back(vec3f(fx, fy, fz));
          }
#endif
        *outOctLowV = oVertice;
        *outOctW    = oWidth;
      }

      void TestOctant::buildOctant(ospray::AMRVolume *amr)
      {
        if (amr == NULL || amr->accel->leaf.size() <= 0)
          return;
        size_t nLeaf = amr->accel->leaf.size();
        printf("Start to generate Octant from %ld leaves!\n", nLeaf);
        std::vector<vec3f> *octLowV = new std::vector<vec3f>[nLeaf];
        std::vector<float> *octW    = new std::vector<float>[nLeaf];

        tasking::parallel_for(nLeaf, [&](int leafID) {
          buildOctantByLeaf(amr->accel,
                            amr->accel->leaf[leafID],
                            &octLowV[leafID],
                            &octW[leafID]);
        });

        std::vector<size_t> leafBegin;
        octNum = 0;
        for (size_t i = 0; i < nLeaf; i++) {
          leafBegin.push_back(octNum);
          octNum += octLowV[i].size();
        }
        octVertices.resize(octNum);
        octWidth.resize(octNum);
        printf("Done generating %ld Octants!\n", octNum);

        tasking::parallel_for(nLeaf, [&](const int leafID) {
          std::copy(octLowV[leafID].begin(),
                    octLowV[leafID].end(),
                    &octVertices[leafBegin[leafID]]);
          std::copy(octW[leafID].begin(),
                    octW[leafID].end(),
                    &octWidth[leafBegin[leafID]]);
        });
        printf("Done Merging %ld Octants!\n", octNum);
      }

      void TestOctant::initOctantValue(ospray::AMRVolume *amr)
      {
        if (this->octNum <= 0)
          return;

        vec3f *octantSampleVetex = (vec3f *)this->octVertices.data();
        float *octWidth          = (float *)this->octWidth.data();
        float *octantValue = new float[this->octVertices.size() * 8];

        std::cout << "Start to reconstruct Octant's value with Octant Method..."
                  << std::endl;
        const size_t blockSize = 1024 * 16;
        size_t blockNum        = this->octVertices.size() / blockSize;
        size_t lastBlock       = this->octVertices.size() % blockSize;

        auto methodStringFromEnv =
            ospcommon::utility::getEnvVar<std::string>("IMPI_AMR_METHOD");

        std::string reconstMethod = methodStringFromEnv.value_or("octant");

        PRINT(reconstMethod);

        tasking::parallel_for(blockNum, [&](int blockID) {
          const size_t begin = size_t(blockID) * blockSize;
          const size_t end =
              std::min(begin + blockSize, this->octVertices.size());
          if (reconstMethod == "octant")
            ispc::getOctantValue_Octant(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
          else if (reconstMethod == "nearst")
            ispc::getOctantValue_Nearst(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
          else if (reconstMethod == "current")
            ispc::getOctantValue_Current(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
        });

        tasking::parallel_for(lastBlock, [&](int taskid) {
          const size_t begin = blockNum * blockSize + taskid;
          if (reconstMethod == "octant")
            ispc::getOctantValue_Octant(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
          else if (reconstMethod == "nearst")
            ispc::getOctantValue_Nearst(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
          else if (reconstMethod == "current")
            ispc::getOctantValue_Current(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
        });

        std::cout << "Get the Value From ISPC..." << std::endl;
        this->octValueBuffer = octantValue;

        std::cout << "Octant's values are set...min:"
                  << *std::min_element(
                         octantValue,
                         octantValue + this->octVertices.size() * 8)
                  << ", max:"
                  << *std::max_element(
                         octantValue,
                         octantValue + this->octVertices.size() * 8)
                  << std::endl;
      }

      TestOctant::~TestOctant()
      {
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
          auto box = box3fa(this->octVertices[i], this->octVertices[i] + vec3f(this->octWidth[i]));
          if (octRange[i].contains(isoValue) && isInClapBox(box)) {
            activeVoxels.push_back(i);
          }
        }
   
        std::cout<<"IsoValue = "<<isoValue<<", ActiveVoxels ="<<activeVoxels.size()<<std::endl;
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        return box3fa(this->octVertices[voxelRef],
                              this->octVertices[voxelRef] + vec3f(this->octWidth[voxelRef]));
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        size_t startIdx = voxelRef * 8;
        voxel.bounds = box3fa(this->octVertices[voxelRef],
                              this->octVertices[voxelRef] + vec3f(this->octWidth[voxelRef]));
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
