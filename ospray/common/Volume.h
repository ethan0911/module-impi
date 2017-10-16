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
#include "ospcommon/range.h"

namespace ospray {
  namespace impi {

    using namespace ospcommon;

    typedef ospcommon::range_t<float> Range;
    
    struct CellRef {
      CellRef() {};
      CellRef(const vec3i &idx)
        : x(idx.x), y(idx.y), z(idx.z)
      {}
      
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

      inline Range getRange() const {
        Range range;
        array3D::for_each(vec3i(2),[&](const vec3i idx) {
            range.extend(vtx[idx.z][idx.y][idx.x]);
          });
        return range;
      }
    };
    
    /*! defines basic logcial abstraction for a volume class we can
        build a iso-dat structure over. should _probably_ at some time
        get merged into some more ospray/common like volume thingy */
    struct LogicalVolume {
      /*! return dimensions (in voxels, not cells!) of the underlying volume */
      virtual vec3i getDims() const = 0;
      
      /* query get the eight corner voxels (in float) of given cell */
      virtual void getCell(Cell &cell, const vec3i &cellIdx) const = 0;

      /* query get the eight corner voxels (in float) of given cell */
      inline Cell getCell(const vec3i &idx) const { Cell c; getCell(c,idx); return c; }
      
      /*! create a list of all the cell references in [lower,upper)
          whose value range overlaps the given iso-value */ 
      virtual void filterVoxelsThatOverLapIsoValue(std::vector<CellRef> &out,
                                                   const vec3i &lower,
                                                   const vec3i &upper,
                                                   const float iso) const = 0;

      /*! build list of all cells that fulfill the given filter lambda */
      template<typename Lambda>
      void filterVoxels(std::vector<CellRef> &out, Lambda filter);
      
      /*! create a list of *all* the cell references in the entire volume
        whose value range overlaps the given iso-value */ 
      void filterAllVoxelsThatOverLapIsoValue(std::vector<CellRef> &out,
                                              const float iso) const;
    };

    /*! load a test data set */
    std::shared_ptr<LogicalVolume> loadTestDataSet();

    /*! defines an _actual_ implementation of a volume */
    template<typename T>
    struct VolumeT : public LogicalVolume, array3D::ActualArray3D<T> {

      /*! constructor */
      VolumeT(const vec3i &dims) : array3D::ActualArray3D<T>(dims) {}
      
      /*! return dimensions (in voxels, not cells!) of the underlying volume */
      virtual vec3i getDims() const override { return this->size(); }

      static std::shared_ptr<LogicalVolume> loadRAW(const std::string fileName,
                                                    const vec3i &dims);
      
      
      inline Range getRangeOfCell(const vec3i &cellIdx) const
      {
        Cell c;
        getCell(c,cellIdx);
        return c.getRange();
      }
      
      /* query get the eight corner voxels (in float) of given cell */
      virtual void getCell(Cell &cell, const vec3i &cellIdx) const override
      {
        array3D::for_each(vec3i(2),[&](const vec3i vtxIdx){
            cell.vtx[vtxIdx.z][vtxIdx.y][vtxIdx.x] = this->get(cellIdx+vtxIdx);
          });
      }
      
      /*! create a list of all the cell references in [lower,upper)
          whose value range overlaps the given iso-value */ 
      virtual void filterVoxelsThatOverLapIsoValue(std::vector<CellRef> &out,
                                                   const vec3i &lower,
                                                   const vec3i &upper,
                                                   const float iso) const override
      {
        array3D::for_each(lower,upper,[&](const vec3i &idx) {
            if (this->getRangeOfCell(idx).contains(iso))
              out.push_back(CellRef(idx));
          });
      }
    };

    /*! build list of all cells that fulfill the given filter lambda */
    template<typename Lambda>
    inline void LogicalVolume::filterVoxels(std::vector<CellRef> &out, Lambda filter)
    {
      /*! for now, do this single-threaded: \todo use tasksys ... */
      out.clear();
      array3D::for_each(getDims(),[&](const vec3i &idx) {
          if (filter(this,idx))
            out.push_back(CellRef(idx));
        });
    }

    
  } // ::ospray::impi
} // ::ospray
