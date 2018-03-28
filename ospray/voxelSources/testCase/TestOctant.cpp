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
#include "ospray/common/Data.h"
#include "compute_voxels_ispc.h"

#include <time.h>
#include <numeric>

#ifndef speedtest__
#define speedtest__(data)                                           \
  for (long blockTime = NULL;					    \
       (blockTime == NULL ? (blockTime = clock()) != NULL : false); \
       printf(data" Calculation Time: %.9fs \n",		    \
	      (double)(clock() - blockTime) / CLOCKS_PER_SEC))
#endif
  
namespace ospray {
  namespace impi {
    namespace testCase {

      /*! constructors and distroctors */
      TestOctant::TestOctant() {}
      TestOctant::~TestOctant()
      {
        if (octValueBuffer != NULL) { delete[] octValueBuffer; }
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        return box3fa(this->octVertices[voxelRef],
		      this->octVertices[voxelRef] + 
		      vec3f(this->octWidth[voxelRef]));
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        size_t startIdx = voxelRef * 8;
        voxel.bounds = box3fa(this->octVertices[voxelRef],
                              this->octVertices[voxelRef] + 
			      vec3f(this->octWidth[voxelRef]));
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

      /*! create lits of *all* voxel (refs) we want to be considered for
       * interesction */
      void TestOctant::getActiveVoxels(std::vector<VoxelRef> &activeVoxels, 
				       float isoValue) const
      {
	// 		
	// compute number of octants first
	//
	const std::unique_ptr<ospray::amr::AMRAccel> &accel = 
	  amrVolumePtr->accel;
        const size_t nLeaf = accel->leaf.size();
        std::vector<size_t> numOfLeafOctants(nLeaf, size_t(0));
	speedtest__("#osp:impi: Compute Number of Octants") {
	  tasking::parallel_for(nLeaf, [&] (size_t lid) {
	      const ospray::amr::AMRAccel::Leaf &lf = accel->leaf[lid];
	      const float fcw = lf.brickList[0]->cellWidth; // cell width
	      const float hcw = 0.5f * fcw;            // half cell width
	      const vec3f &lower = lf.bounds.lower;
	      const vec3f &upper = lf.bounds.upper;
	      const size_t nx = std::round((upper.x - lower.x) / hcw);
	      const size_t ny = std::round((upper.y - lower.y) / hcw);
	      const size_t nz = std::round((upper.z - lower.z) / hcw);
	      // add inner cells
	      numOfLeafOctants[lid] += 
		(nx / size_t(2) - size_t(1)) * 
		(ny / size_t(2) - size_t(1)) * 
		(nz / size_t(2) - size_t(1));
	      // bottom top boundray cells
	      numOfLeafOctants[lid] += ny * nx * 2;
	      // left right boundray cells
	      numOfLeafOctants[lid] += nz * ny * 2;
	      // front back boundary cells
	      numOfLeafOctants[lid] += nz * nx * 2;
	    });
	}
	const size_t numOfOctants = 
	  std::accumulate(numOfLeafOctants.begin(),
			  numOfLeafOctants.end(),
			  size_t(0));
	printf("#osp:impi: numOfOctants %zu\n", numOfOctants);
        if (numOfOctants <= 0) {
	  throw std::runtime_error("No octants are found");
	}
	//
	// safty check
	//
	assert(numOfOctants == octNum);
	//
	// iterate through all voxels
	//
        const size_t blockSize = 1024 * 16;
        const size_t blockNum  = this->octVertices.size() / blockSize;
        const size_t lastBlock = this->octVertices.size() % blockSize;
        const auto methodStringFromEnv =
            ospcommon::utility::getEnvVar<std::string>("IMPI_AMR_METHOD");
        const std::string method = methodStringFromEnv.value_or("octant");
	printf("#osp:impi: method %s\n", method.c_str());
	//
	//
	//
#if 0
        tasking::parallel_for(blockNum, [&](int blockID) {
	    const size_t begin = size_t(blockID) * blockSize;
	      ispc::getOctantValue_Octant(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                blockSize);
        });

        tasking::parallel_for(lastBlock, [&](int taskid) {
          const size_t begin = blockNum * blockSize + taskid;
	  ispc::getOctantValue_Octant(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
        });

	// ------
#endif
        activeVoxels.clear(); // the output
        for (size_t i = 0; i < this->octNum; i++) {
          auto box = box3fa(this->octVertices[i], 
			    this->octVertices[i] + 
			    vec3f(this->octWidth[i]));
          if (octRange[i].contains(isoValue) && isInClapBox(box)) {
            activeVoxels.push_back(i);
          }
        }
   
        std::cout << "#osp:impi: " 
		  << "IsoValue = " << isoValue << " "
		  << "Number of ActiveVoxels = "
		  << activeVoxels.size() 
		  << std::endl;
	
      }

      void TestOctant::initOctant(ospray::AMRVolume * amr)
      {
	//------------------------------------------------------------------//
	// check if voxels are valid
	//------------------------------------------------------------------//
        if (!amr) {
	  throw std::runtime_error("Empty amr volume");	  
	}
        if (amr->accel->leaf.size() <= 0) {
          throw std::runtime_error("AMR Volume has no leaf");
	}
	amrVolumePtr = amr;

	//return; // Qi: hijack here
	
	//------------------------------------------------------------------//
#if 1	
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
        // b1.upper.z *= 0.5f;
        this->clapBoxes.push_back(b1);

        speedtest__("Speed: ")
        {
          this->octRange.resize(this->octNum);
          for (size_t i = 0; i < this->octNum; i++) {
            getOctrange(i, &this->octRange[i]);
          }
        }
        std::cout << "Done Init Octant Value!" << std::endl;
#endif
      }

      void 
      TestOctant::buildOctantByLeaf(std::unique_ptr<ospray::amr::AMRAccel> 
				    &accel,
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
        *outOctLowV = oVertice; // left, down, vertex
        *outOctW    = oWidth;   // width
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

	speedtest__("Speed Compute Values: ") {
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
          else if (reconstMethod == "nearest")
            ispc::getOctantValue_Nearest(
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
          else if (reconstMethod == "nearest")
            ispc::getOctantValue_Nearest(amr->getIE(),
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
	}

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



    }
  }
}
