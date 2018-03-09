// ======================================================================== //
// Copyright SCI Institute, University of Utah, 2018
// ======================================================================== //

#pragma once

#include "ospray/ospray.h"
#include "ospcommon/malloc.h"
#include "ospcommon/math.h"
#include "ospcommon/vec.h"
#include "ospcommon/box.h"
#include "ospcommon/range.h"
#include "common/xml/XML.h"
#include <vector>
#include <string>

namespace ospray {

  namespace amr {

    //! AMR SG node with Chombo style structure
    struct AMRVolume
    {
      AMRVolume() : maxLevel(1 << 30), amrMethod("current") {}
      ~AMRVolume() {
	for (auto *ptr : brickPtrs) {
	  delete [] ptr;
	}
	for (auto& obj : brickData) {
	  ospRelease(obj);
	}
	if (volume != nullptr) {
	  ospRelease(volume);
	  ospRelease(ospBrickData);
	  ospRelease(ospBrickInfo);
	}
      };

      void Load(const xml::Node &node);
      OSPVolume Create(OSPTransferFunction tfn) {

	volume = ospNewVolume("amr_volume");

	for (size_t bID = 0; bID < brickInfo.size(); bID++) {
	  const auto &bi = brickInfo[bID];
	  OSPData data = ospNewData(bi.size().product(),
				    OSP_FLOAT,
				    this->brickPtrs[bID],
				    OSP_DATA_SHARED_BUFFER);
	  this->brickData.push_back(data);
	}

	ospBrickData = ospNewData(brickData.size(),
				  OSP_DATA,
				  (OSPObject*)(brickData.data()),
				  OSP_DATA_SHARED_BUFFER);
	
	ospBrickInfo = ospNewData(brickInfo.size() * sizeof(BrickInfo),
				  OSP_UCHAR,
				  (uint8_t*)(brickInfo.data()),
				  OSP_DATA_SHARED_BUFFER);
	
	ospSetData(volume, "brickData", ospBrickData);
	ospSetData(volume, "brickInfo", ospBrickInfo);
	ospSetString(volume, "voxelType", "float");
	
	ospSetObject(volume, "transferFunction", tfn);
	ospSetVec2f(volume, "voxelRange", (const osp::vec2f&)voxelRange);
	ospSet1i(volume, "gradientShadingEnabled", 0);
	ospSet1i(volume, "preIntegration", 0);
	ospSet1i(volume, "singleShade", 0);
	ospSet1i(volume, "adaptiveSampling", 0);
	ospSet1f(volume, "samplingRate", 1.f);

	ospCommit(volume);
	return volume;

      }

      // ------------------------------------------------------------------
      // this is the way we're passing over the data. for each input
      // box we create one data array (for the data values), and one
      // 'BrickInfo' descriptor that describes it. we then pass two
      // arrays - one array of all AMRBox'es, and one for all data
      // object handles. order of the brickinfos MUST correspond to
      // the order of the brickDataArrays; and the brickInfos MUST be
      // ordered from lowest level to highest level (inside each level
      // it doesn't matter)
      // ------------------------------------------------------------------
      struct BrickInfo
      {
        ospcommon::box3i box;
        int level;
        float dt;

        ospcommon::vec3i size() const
        {
          return box.size() + ospcommon::vec3i(1);
        }
      };

      // OSPRay data handlers
      OSPVolume volume = nullptr;
      OSPData ospBrickData = nullptr;
      OSPData ospBrickInfo = nullptr;

      // ID of the data component we want to render (each brick can
      // contain multiple components)
      int componentID{0};
      int maxLevel;
      ospcommon::vec2f voxelRange; // the copy of value range, which will be updated in Load
      ospcommon::range1f valueRange;
      ospcommon::box3f bounds;
      std::string amrMethod;
      std::vector<OSPData> brickData;
      std::vector<BrickInfo> brickInfo;
      std::vector<float *> brickPtrs;

    };

  };

  namespace ParseOSP {
    std::shared_ptr<ospray::amr::AMRVolume> loadOSP(const std::string &fileName);
  };
  
};
