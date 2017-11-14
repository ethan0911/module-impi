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

namespace ospray {
  namespace impi { 
    namespace testCase {

      TestOctant::TestOctant(){
        //std::string fileName = "Octant1.oct";
        //parseOctant(fileName);
        //PRINT(this->octants.size());
      }

      // void TestOctant::parseOctant(std::string fileName){
      //   std::ifstream inFile;
      //   inFile.open(fileName);
      //   if(!inFile){
      //     return;
      //     //perror(("Error while opening file" + fileName).c_str());
      //   }
  
      //   std::string sOctantNum;
      //   getline(inFile,sOctantNum);
      //   int octantNum = stoi(sOctantNum);

      //   std::string line;
      //   while(std::getline(inFile,line)){
      //     if(line == "")
      //       break;
      //     std::stringstream ss(line);
      //     vec3f p;
      //     float cWidth;
      //     float value[8];
      //     ss >> p.x >> p.y >> p.z >> cWidth >> value[0] >> value[1] >>
      //         value[2] >> value[3] >> value[4] >> value[5] >> value[6] >>
      //         value[7];
      //     testCase::Octant oct;
      //     oct.bounds.lower = p;
      //     oct.width = cWidth;
      //     oct.bounds.upper = p + vec3f(oct.width);
      //     oct.vertexValue[0][0][0]    = value[0];
      //     oct.vertexValue[0][0][1]    = value[1];
      //     oct.vertexValue[0][1][0]    = value[2];
      //     oct.vertexValue[0][1][1]    = value[3];
      //     oct.vertexValue[1][0][0]    = value[4];
      //     oct.vertexValue[1][0][1]    = value[5];
      //     oct.vertexValue[1][1][0]    = value[6];
      //     oct.vertexValue[1][1][1]    = value[7];
      //     octants.push_back(oct);
      //   }
      //   if(inFile.bad())
      //     perror(("Error while reading file" + fileName).c_str());
      //   inFile.close();

      //   // PRINT(octantNum);
      // }
      void TestOctant::initData(int octNum,
                                vec3f *octLowerPnt,
                                float *octWidth,
                                float *octantValue)
      {
        for (int i = 0; i < octNum; ++i) {
          testCase::Octant oct;
          oct.bounds.lower         = octLowerPnt[i];
          oct.width                = octWidth[i];
          oct.bounds.upper         =  oct.bounds.lower  + vec3f(oct.width);
          int index = 8*i;
          oct.vertexValue[0][0][0] = octantValue[index++];
          oct.vertexValue[0][0][1] = octantValue[index++];
          oct.vertexValue[0][1][0] = octantValue[index++];
          oct.vertexValue[0][1][1] = octantValue[index++];
          oct.vertexValue[1][0][0] = octantValue[index++];
          oct.vertexValue[1][0][1] = octantValue[index++];
          oct.vertexValue[1][1][0] = octantValue[index++];
          oct.vertexValue[1][1][1] = octantValue[index++];
          octants.push_back(oct);
        }
      }

      /*! create lits of *all* voxel (refs) we want to be considered for
       * interesction */
      void TestOctant::getActiveVoxels(std::vector<VoxelRef> &activeVoxels, float isoValue) const 
      {
        activeVoxels.clear();
        for (int i=0;i<octants.size();i++)
          activeVoxels.push_back(i);
      }

      /*! compute world-space bounds for given voxel */
      box3fa TestOctant::getVoxelBounds(const VoxelRef voxelRef) const 
      {
        const Octant& oct = octants[(const int)voxelRef];
        return box3fa(oct.bounds.lower,oct.bounds.upper);
      }

      /*! get full voxel - bounds and vertex values - for given voxel */
      Impi::Voxel TestOctant::getVoxel(const VoxelRef voxelRef) const 
      {
        Impi::Voxel voxel;
        const Octant& oct = octants[(const int)voxelRef];
        voxel.bounds = box3fa(oct.bounds.lower,oct.bounds.upper);
        array3D::for_each(vec3i(2),[&](const vec3i vtx){
          voxel.vtx[vtx.z][vtx.y][vtx.x] = oct.vertexValue[vtx.z][vtx.y][vtx.x];
        });
        return voxel;
      }
    }
  }
}
