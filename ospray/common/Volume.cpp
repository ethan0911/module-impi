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

#include "Volume.h"

namespace ospray {
  namespace impi {

    /*! create a list of *all* the cell references in the entire volume
      whose value range overlaps the given iso-value 
    */ 
    void LogicalVolume::filterAllVoxelsThatOverLapIsoValue(std::vector<CellRef> &out,
                                              const float iso) const
    {
      /*! for now, do this single-threaded: \todo use tasksys ... */
      out.clear();
      filterVoxelsThatOverLapIsoValue(out,vec3i(0),getDims()-vec3i(1),iso);
    }


    /* geneate a simple linear basis function centerred around cente,r
       with max radius r, and center wieght 1, and evalute it for
       'pos' */
    inline float genBlob(const vec3f &pos, const vec3f &center, const float r)
    {
      float dist = length(pos-center);
      if (dist > r) return 0.f;
      return 1.f - (dist/r);
    }
    
    std::shared_ptr<LogicalVolume> loadTestDataSet()
    {
#if 1
      const vec3i dims(64);
      std::shared_ptr<VolumeT<float>> vol = std::make_shared<VolumeT<float>>(dims);

      array3D::for_each(dims,[&](const vec3i &idx){
          float val = 0.f;
          vec3f pos = (vec3f(idx)+vec3f(.5f)) * rcp(vec3f(dims));
          val += genBlob(pos,vec3f(.2,.1,.7),.4);
          val += genBlob(pos,vec3f(.3,.3,.2),.3);
          val += genBlob(pos,vec3f(.8,.4,.9),.2);
          val += genBlob(pos,vec3f(.5,.5,.5),.5);
          vol->set(idx,val);
        });
      return vol;
#else
      const std::string fileName = "~/models/magnetic-512-volume/magnetic-512-volume.raw";
      const vec3i dims(512);
      
      FILE *file = fopen(fileName.c_str(),"rb");
      if (!file)
        throw std::runtime_error("could not load volume '"+fileName+"'");
      std::shared_ptr<VolumeT<float>> vol = std::make_shared<VolumeT<float>>(dims);
      size_t numRead = fread(vol->value,dims.product(),sizeof(float),file);
      assert(numRead == dims.product());
      return vol;
#endif
    }

  } // ::ospray::impi
} // ::ospray

    
