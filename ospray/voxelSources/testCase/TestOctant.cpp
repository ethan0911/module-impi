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
#include "common/Data.h"
#include "compute_voxels_ispc.h"

#include <time.h>
#include <numeric>

#include <mutex>
std::mutex lock;

#ifndef speedtest__
#define speedtest__(data)                                           \
  for (long blockTime = NULL;					    \
       (blockTime == NULL ? (blockTime = clock()) != NULL : false); \
       printf(data" Calculation Time: %.9fs \n",		    \
	      (double)(clock() - blockTime) / CLOCKS_PER_SEC))
#endif

#define ON_THE_FLY 1

namespace ospray {
  namespace impi {
    namespace testCase {

      // extern "C" void  externC_push_back(void * _c_vector,
      // 					 void * _c_ptr,
      // 					 const float v0,
      // 					 const float v1,
      // 					 const float v2,
      // 					 const float v3,
      // 					 const float v4,
      // 					 const float v5,
      // 					 const float v6,
      // 					 const float v7,
      // 					 const float c0,
      // 					 const float c1,
      // 					 const float c2,
      // 					 const float cellwidth)
      // {	
      // 	auto c_ptr = (TestOctant*) _c_ptr;
      // 	vec3f coordinate(c0, c1, c2);
      // 	if (c_ptr->inClipBox(box3fa(coordinate, coordinate + cellwidth))) {  
      // 	  auto c_vector = (std::vector<Octant>*) _c_vector;
      // 	  c_vector->emplace_back();	
      // 	  c_vector->back().vtx[0][0][0] = v0;
      // 	  c_vector->back().vtx[0][0][1] = v1;
      // 	  c_vector->back().vtx[0][1][0] = v2;
      // 	  c_vector->back().vtx[0][1][1] = v3;
      // 	  c_vector->back().vtx[1][0][0] = v4;
      // 	  c_vector->back().vtx[1][0][1] = v5;
      // 	  c_vector->back().vtx[1][1][0] = v6;
      // 	  c_vector->back().vtx[1][1][1] = v7;
      // 	  c_vector->back().lowerleft = coordinate;
      // 	  c_vector->back().cellwidth = cellwidth;
      // 	}
      // }
      extern "C" void  externC_push_back(void * _c_vector,
					 void * _c_ptr,
					 const uint32_t lid,
					 const uint32_t oid,
					 const float c0,
					 const float c1,
					 const float c2,
					 const float cellwidth)
      {	
	const auto c_ptr = (TestOctant*) _c_ptr;
	const auto box = box3fa(vec3f(c0, c1, c2), 
				vec3f(c0, c1, c2) + cellwidth);
	if (c_ptr->inClipBox(box)) {
	  auto c_vector = (std::vector<uint64_t>*) _c_vector;
	  uint64_t idx = (uint64_t(lid) << 32) | uint64_t(oid);
	  c_vector->push_back(idx);
	}
      }


      /*! constructors and distroctors */
      TestOctant::TestOctant() {}
      TestOctant::~TestOctant()
      {
        if (octValueBuffer != NULL) { delete[] octValueBuffer; }
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        // return box3fa(octants[voxelRef].lowerleft,
	// 	      octants[voxelRef].lowerleft + 
	// 	      octants[voxelRef].cellwidth);
	return getVoxel(voxelRef).bounds;
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        // voxel.bounds = box3fa(octants[voxelRef].lowerleft,
	// 		      octants[voxelRef].lowerleft + 
	// 		      octants[voxelRef].cellwidth);        
	// voxel.vtx[0][0][0] = octants[voxelRef].vtx[0][0][0];
        // voxel.vtx[0][0][1] = octants[voxelRef].vtx[0][0][1];
        // voxel.vtx[0][1][0] = octants[voxelRef].vtx[0][1][0];
        // voxel.vtx[0][1][1] = octants[voxelRef].vtx[0][1][1];
        // voxel.vtx[1][0][0] = octants[voxelRef].vtx[1][0][0];
        // voxel.vtx[1][0][1] = octants[voxelRef].vtx[1][0][1];
        // voxel.vtx[1][1][0] = octants[voxelRef].vtx[1][1][0];
        // voxel.vtx[1][1][1] = octants[voxelRef].vtx[1][1][1];
	union {
	  uint32_t out[2];
	  uint64_t in;
	} caster;
	caster.in = voxelRef;
	uint32_t oid = caster.out[0];
	uint32_t lid = caster.out[1];
	// std::cout << lid << " " << oid << std::endl;; 

	const ospray::amr::AMRAccel::Leaf &lf = amrVolumePtr->accel->leaf[lid];	
	const float w = lf.brickList[0]->cellWidth; // cell width
	const float s = lf.brickList[0]->gridToWorldScale;
	const vec3f &lower = lf.bounds.lower;
	const vec3f &upper = lf.bounds.upper;
	const size_t nx = std::round((upper.x - lower.x) * s);
	const size_t ny = std::round((upper.y - lower.y) * s);
	const size_t nz = std::round((upper.z - lower.z) * s);
	// add inner cells
	const auto n1 = 
	  (nx - size_t(1)) * (ny - size_t(1)) * (nz - size_t(1));
	// bottom top boundray cells
	const auto n2 = size_t(8) * ny * nx;
	// left right boundray cells
	const auto n3 = size_t(8) * nz * ny;
	// front back boundary cells
	const auto n4 = size_t(8) * nz * nx;
	//
	float cellwidth;
	ispc::getVoxel(amrVolumePtr->getIE(),			      
		       &voxel.vtx[0][0][0],
		       (ispc::vec3f&)voxel.bounds.lower,
		       cellwidth,
		       w, oid, 
		       (ispc::vec3f&)lower, (ispc::vec3f&)upper,
		       (uint32_t)nx, (uint32_t)ny, (uint32_t)nz,
		       (uint32_t)n1,
		       (uint32_t)(n2 + n1),
		       (uint32_t)(n3 + n2 + n1));
	voxel.bounds.upper = voxel.bounds.lower + cellwidth;
        return voxel;
      }

      /*! create lits of *all* voxel (refs) we want to be considered for
       * interesction */
      void TestOctant::buildActiveVoxels(std::vector<VoxelRef> &activeVoxels, 
					 const float isoValue)
      {
	// octants.clear();
	activeVoxels.clear(); // the output
#if ON_THE_FLY /* Qi: hijack here for on-the-fly calculation */
	//------------------------------------------------------------------//
	// 		
	// initialization
	//
	const auto &accel = amrVolumePtr->accel;
        const auto nLeaf = accel->leaf.size();
        const auto methodFromEnv =
            ospcommon::utility::getEnvVar<std::string>("IMPI_AMR_METHOD");
        const std::string method = methodFromEnv.value_or("octant");
	printf("#osp:impi: method %s\n", method.c_str());

	// //
	// // compute number of octants first
	// //
	// auto *numOfLeafOctants = new size_t[nLeaf]();
	// speedtest__("#osp:impi: Compute Octants") {	  
	//   tasking::parallel_for(nLeaf, [&] (size_t lid) {
	//       //
	//       // meta data
	//       //
	//       const ospray::amr::AMRAccel::Leaf &lf = accel->leaf[lid];
	//       const float w = lf.brickList[0]->cellWidth; // cell width
	//       const float s = lf.brickList[0]->gridToWorldScale;
	//       const vec3f &lower = lf.bounds.lower;
	//       const vec3f &upper = lf.bounds.upper;
	//       const size_t nx = std::round((upper.x - lower.x) * s);
	//       const size_t ny = std::round((upper.y - lower.y) * s);
	//       const size_t nz = std::round((upper.z - lower.z) * s);
	//       //
	//       // number of octants
	//       //
	//       // add inner cells
	//       const auto n1 = 
	// 	(nx - size_t(1)) * (ny - size_t(1)) * (nz - size_t(1));
	//       // bottom top boundray cells
	//       const auto n2 = size_t(8) * ny * nx;
	//       // left right boundray cells
	//       const auto n3 = size_t(8) * nz * ny;
	//       // front back boundary cells
	//       const auto n4 = size_t(8) * nz * nx;
	//       // total
	//       auto &N = numOfLeafOctants[lid];	      
	//       N += n1 + n2 + n3 + n4;
	//     });
	// }
	// size_t numOfOctants(0);
	// for (int i = 0; i < nLeaf; ++i) {
	//   const auto tmp = numOfLeafOctants[i];
	//   numOfLeafOctants[i] = numOfOctants;
	//   numOfOctants += tmp;
	// }
	// //
	// // safty check
	// //
	// printf("#osp:impi: numOfOctants %zu\n", numOfOctants);
        // if (numOfOctants <= 0) {
	//   throw std::runtime_error("No octants are found");
	// }

	//
	// Testing my implementation
	//
        std::cout << "#osp:impi: Number of Leaves " << nLeaf << std::endl;
        std::cout << "#osp:impi: Computing Values Values " << std::endl;
	//auto leafActiveOctants = new std::vector<Octant>[nLeaf];
	auto leafActiveOctants = new std::vector<uint64_t>[nLeaf];
	// repeat the computation here for testing purpose
	speedtest__("#osp:impi: Compute Octants Values") {	  
	  tasking::parallel_for(nLeaf, [&] (size_t lid) {
	      //
	      // meta data
	      //
	      const ospray::amr::AMRAccel::Leaf &lf = accel->leaf[lid];
	      const float w = lf.brickList[0]->cellWidth; // cell width
	      const float s = lf.brickList[0]->gridToWorldScale;
	      const vec3f &lower = lf.bounds.lower;
	      const vec3f &upper = lf.bounds.upper;
	      //TODO: this is wrong, why ??
	      const size_t nx = std::round((upper.x - lower.x) * s);
	      const size_t ny = std::round((upper.y - lower.y) * s);
	      const size_t nz = std::round((upper.z - lower.z) * s);
	      // {
	      // 	const size_t _nx = lf.brickList[0]->dims.x;
	      // 	const size_t _ny = lf.brickList[0]->dims.y;
	      // 	const size_t _nz = lf.brickList[0]->dims.z;
	      // 	lock.lock();
	      // 	std::cout << "nx: " << nx << " "
	      // 		  << "_nx: " << _nx << std::endl;
	      // 	lock.unlock();
	      // }
	      const auto &rg = lf.valueRange;
	      //
	      // number of octants
	      //
	      // add inner cells
	      const auto n1 = 
		(nx - size_t(1)) * (ny - size_t(1)) * (nz - size_t(1));
	      // bottom top boundray cells
	      const auto n2 = size_t(8) * ny * nx;
	      // left right boundray cells
	      const auto n3 = size_t(8) * nz * ny;
	      // front back boundary cells
	      const auto n4 = size_t(8) * nz * nx;
	      // total number
	      const auto N = (n1+n2+n3+n4);
	      //
	      // check leaf range
	      //
	      // lock.lock();
	      // std::cout << lid << " " 
	      // 		<< &leafActiveOctants[lid] << " " 
	      // 		<< leafActiveOctants[lid].size() << std::endl;
	      // lock.unlock();
	      // const size_t blockSize(32 * 16);
	      // const size_t blockNum = (N + blockSize - 1) / blockSize;
	      // const size_t off = numOfLeafOctants[lid];
	      // tasking::parallel_for(blockNum, [&](size_t blockID) {
	      // const size_t b = blockID * blockSize;
	      // const size_t e = std::min(b + blockSize, N);		 
	      const size_t b = 0;
	      const size_t e = N;
	      ispc::getVoxel_Octant(amrVolumePtr->getIE(),
				    this, &leafActiveOctants[lid],
				    isoValue, w, lid, 
				    (ispc::vec3f&)lower, (ispc::vec3f&)upper,
				    (uint32_t)b,  (uint32_t)e, 
				    (uint32_t)nx, (uint32_t)ny, (uint32_t)nz,
				    (uint32_t)n1,
				    (uint32_t)(n2 + n1),
				    (uint32_t)(n3 + n2 + n1));
	      //});
	    });
	}
	//
	//
	//
	// delete[] numOfLeafOctants;
        std::cout << "#osp:impi: Done Computing Values Values" << std::endl;

	std::vector<size_t> begin(nLeaf, size_t(0));
	size_t n(0);
	for (int lid = 0; lid < nLeaf; ++lid) {
	  begin[lid] = n;
	  n += leafActiveOctants[lid].size();
	  // std::cout << leafActiveOctants[lid].size() << " "
	  // 	    << begin[lid] << std::endl;
	}
	// octants.resize(n);
	activeVoxels.resize(n);
        tasking::parallel_for(nLeaf, [&](const size_t lid) {
          std::copy(leafActiveOctants[lid].begin(),
                    leafActiveOctants[lid].end(),
                    &activeVoxels[begin[lid]]);
        });

	delete[] leafActiveOctants;

        // std::cout << "#osp:impi: Hack: Octant's values are set... min: "
        //           << *std::min_element(octVBuff,
	// 			       octVBuff +
	// 			       numOfOctants * 8)
        //           << ", max: "
        //           << *std::max_element(octVBuff,
	// 			       octVBuff + 
	// 			       numOfOctants * 8)
        //           << std::endl;
	
        std::cout << "Done Init Octant Value! " << activeVoxels.size() << std::endl;
	//------------------------------------------------------------------//	
#endif
        //activeVoxels.clear(); // the output
	//for (int i = 0; i < octants.size(); ++i)  {
	//  activeVoxels.push_back(i);
	//}
      }

      void TestOctant::getActiveVoxels(std::vector<VoxelRef> &activeVoxels, 
				       float isoValue) const
      {
	//
	// clean
	//
        activeVoxels.clear(); // the output

        for (size_t i = 0; i < this->octNum; i++) {
          auto box = box3fa(this->octVertices[i], 
			    this->octVertices[i] + 
			    vec3f(this->octWidth[i]));
          if (octRange[i].contains(isoValue) && inClipBox(box)) {
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

	//------------------------------------------------------------------//
	// compute default bbox
	//------------------------------------------------------------------//
        const auto wb = amr->accel->worldBounds;
        clipBoxes.push_back(box3fa(amr->accel->worldBounds.lower, 
				   amr->accel->worldBounds.upper));
	
#if !(ON_THE_FLY)
	//------------------------------------------------------------------//
        std::cout << "Start to Init Octant Value" << std::endl;
        
        buildOctant(amr);
        initOctantValue(amr);

        auto wb = amr->accel->worldBounds;

        // vec3f center = wb.center();
        // box3fa b1 = box3fa(wb.lower,wb.upper);
        // b1.upper.x *= 0.5f;
        // box3fa b2 = box3fa(vec3f(center.x,0.f,0.f),
	// 		      vec3f(wb.upper.x,center.y,center.z));
        // box3fa b3 = box3fa(vec3f(center.x,center.y,0.f),wb.upper);
        // this->clipBoxes.push_back(b1);
        // this->clipBoxes.push_back(b2);
        // this->clipBoxes.push_back(b3);

        box3fa box = box3fa(wb.lower, wb.upper); // box.upper.z *= 0.5f;
        this->clipBoxes.push_back(box);

        speedtest__("Speed: ")
        {
          this->octRange.resize(this->octNum);
          for (size_t i = 0; i < this->octNum; i++) {
            getOctrange(i, &this->octRange[i]);
          }
        }
        std::cout << "Done Init Octant Value!" << std::endl;
	//------------------------------------------------------------------//
#endif
      }

      //---------------------------------------------------------------------------//
      //
      //---------------------------------------------------------------------------//
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
        // bottom top boundray cells
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
        // finest VA
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
	// Qi: what is this ?
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

        std::cout << "Start to reconstruct Octant's value with Octant Method"
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
            ispc::getAMRValue_Octant(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
          else if (reconstMethod == "nearest")
            ispc::getAMRValue_Nearest(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
          else if (reconstMethod == "current")
            ispc::getAMRValue_Current(
                amr->getIE(),
                &octantValue[begin * 8],
                (ispc::vec3f *)&octantSampleVetex[begin],
                &octWidth[begin],
                end - begin);
        });

        tasking::parallel_for(lastBlock, [&](int taskid) {
          const size_t begin = blockNum * blockSize + taskid;
          if (reconstMethod == "octant")
            ispc::getAMRValue_Octant(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
          else if (reconstMethod == "nearest")
            ispc::getAMRValue_Nearest(amr->getIE(),
                                 &octantValue[begin * 8],
                                 (ispc::vec3f *)&octantSampleVetex[begin],
                                 &octWidth[begin],
                                 1);
          else if (reconstMethod == "current")
            ispc::getAMRValue_Current(amr->getIE(),
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
