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
    
    inline std::shared_ptr<LogicalVolume> loadTestDataSet()
    {
      const std::string fileName = "~/models/magnetic-512-volume/magnetic-512-volume.raw";
      const vec3i dims(512);
      
      FILE *file = fopen(fileName.c_str(),"rb");
      if (!file)
        throw std::runtime_error("could not load volume '"+fileName+"'");
      std::shared_ptr<VolumeT<float>> vol = std::make_shared<VolumeT<float>>(dims);
      size_t numRead = fread(vol->value,dims.product(),sizeof(float),file);
      assert(numRead == dims.product());
      return vol;
    }

  } // ::ospray::impi
} // ::ospray

    
