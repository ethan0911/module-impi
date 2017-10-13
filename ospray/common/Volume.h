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

#include "ospcommon/array3D/Array3D.h"
#include "ospcommon/array3D/for_each.h"

namespace ospray {
  namespace impi {

    using namespace ospcommon;

    struct VoxelRef {
      VoxelRef() {};
      VoxelRef(const vec3i &idx) : x(idx.x), y(idx.y), z(idx.z) {}
      
      union {
        uint64_t id;
        struct {
          uint64_t x:21;
          uint64_t y:21;
          uint64_t z:21;
        };
      };
    };

    /*! a cell with 8 corner voxels */
    struct Cell {
      float vtx[2][2][2];
    };
    
    /*! defines basic logcial abstraction for a volume class we can
        build a iso-dat structure over. should _probably_ at some time
        get merged into some more ospray/common like volume thingy */
    struct LogicalVolume {
      /*! return dimensions (in voxels, not cells!) of the underlying volume */
      virtual vec3i getDims() const override = 0;
      
      /* query get the eight corner voxels (in float) of given cell */
      virtual void getCell(Cell &cell, const vec3i &cellIdx);
      
      /*! create a list of all the cell references in [lower,upper)
          whose value range overlaps the given iso-value */ 
      virtual void filterVoxelsThatOverLapIsoValue(std::vector<VoxelRef> &out,
                                                   const vec3i &lower,
                                                   const vec3i &upper,
                                                   const float iso) const override = 0;
    };

    /*! defines an _actual_ implementation of a volume */
    template<typename T>
    struct VolumeT : public AbstractVolume, public array3D::ActualArray3D<T> {
      /*! return dimensions (in voxels, not cells!) of the underlying volume */
      virtual vec3i getDims() const override ;
      
      /* query get the eight corner voxels (in float) of given cell */
      virtual void getCell(Cell &cell, const vec3i &cellIdx);
      {
        array3D::for_each(vec3i(2),[&](const vec3i vtxIdx){
            cell.v[vtxIdx.z][vtxIdx.y][vtxIdx.x] = this->get(idx+vtxIdx);
          });
      }
      
      /*! create a list of all the cell references in [lower,upper)
          whose value range overlaps the given iso-value */ 
      virtual void filterVoxelsThatOverLapIsoValue(std::vector<VoxelRef> &out,
                                                   const vec3i &lower,
                                                   const vec3i &upper,
                                                   const float iso) const override
      {
        array3D::for_each(lower,upper,[&](const vec3i &idx) {
            if (this->getRange(idx).contains(iso))
              out.push_back(VoxelRef(idx));
          });
      }
    };

  }
}
