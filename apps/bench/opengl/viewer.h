#pragma once
#ifndef OSPRAY_VIEWER_H
#define OSPRAY_VIEWER_H

#include "ospray/ospray.h"

namespace ospray {
  namespace viewer {
    int Init(const int ac, const char** av);
    void Render(int);
    void Handler(OSPCamera);
    void Handler(OSPTransferFunction);
    void Handler(OSPRenderer);
  };
};


#endif//OSPRAY_VIEWER_H
