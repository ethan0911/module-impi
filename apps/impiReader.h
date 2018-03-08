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

  namespace impi {

    //! AMR SG node with Chombo style structure
    struct AMRVolume
    {
      AMRVolume() : maxLevel(1 << 30), amrMethod("current") {}
      ~AMRVolume() {
	for (auto *ptr : brickPtrs)
	  delete [] ptr;
      };

      void Load(const xml::Node &node);

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

      // ID of the data component we want to render (each brick can
      // contain multiple components)
      int componentID{0};
      int maxLevel;
      ospcommon::range1f valueRange;
      std::vector<OSPData> brickData;
      std::vector<BrickInfo> brickInfo;
      std::vector<float *> brickPtrs;
      ospcommon::box3f bounds;
      std::string amrMethod;

    };

  };

};
