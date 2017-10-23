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

#include <vector>

#include "common/sg/SceneGraph.h"
#include "common/sg/Renderer.h"
#include "common/sg/common/Data.h"
#include "common/sg/geometry/Geometry.h"

#include "CommandLine.h"

#include "exampleViewer/widgets/imguiViewer.h"

/*! _everything_ in the ospray core universe should _always_ be in the
  'ospray' namespace. */
namespace ospray {

  /*! though not required, it is good practice to put any module into
    its own namespace (isnide of ospray:: ). Unlike for the naming of
    library and init function, the naming for this namespace doesn't
    particularlly matter. E.g., 'impi', 'module_blp',
    'bilinar_patch' etc would all work equally well. */
  namespace impi {

    /*! A Simple Triangle Mesh that stores vertex, normal, texcoord,
        and vertex color in separate arrays */
    struct ImpiSGNode : public sg::Geometry
    {
      ImpiSGNode() : Geometry("impi") {}

      box3f bounds() const override
      {
        return box3f(vec3f(0),vec3f(1));
      }
      void setFromXML(const xml::Node &node,
                      const unsigned char *binBasePtr) override {}
    };

    // use ospcommon for vec3f etc
    using namespace ospcommon;

    extern "C" int main(int ac, const char **av)
    {
      int init_error = ospInit(&ac, av);
      if (init_error != OSP_NO_ERROR) {
        std::cerr << "FATAL ERROR DURING INITIALIZATION!" << std::endl;
        return init_error;
      }

      auto device = ospGetCurrentDevice();
      if (device == nullptr) {
        std::cerr << "FATAL ERROR DURING GETTING CURRENT DEVICE!" << std::endl;
        return 1;
      }

      ospDeviceSetStatusFunc(device, [](const char *msg) { std::cout << msg; });
      ospDeviceSetErrorFunc(device,
                            [](OSPError e, const char *msg) {
                              std::cout << "OSPRAY ERROR [" << e << "]: "
                                        << msg << std::endl;
                              std::exit(1);
                            });

      ospDeviceCommit(device);

      // access/load symbols/sg::Nodes dynamically
      loadLibrary("ospray_sg");
      ospLoadModule("impi");

      ospray::imgui3D::init(&ac,av);

      // parse the commandline; complain about anything we do not
      // recognize
      CommandLine args(ac,av);

      auto renderer_ptr = sg::createNode("renderer", "Renderer");
      auto &renderer = *renderer_ptr;

      auto &win_size = ospray::imgui3D::ImGui3DWidget::defaultInitSize;
      renderer["frameBuffer"]["size"] = win_size;

      
      // renderer["rendererType"] = std::string("raycast_Ns");
      renderer["rendererType"] = std::string("scivis");

      auto &world = renderer["world"];

      // auto &patchesInstance = world.createChild("patches", "Instance");

#if 0
      ????
      const std::string fileName = "test.vol";
      auto importerNode_ptr = sg::createNode(ss.str(), "Importer")->nodeAs<sg::Importer>();;
      auto &importerNode = *importerNode_ptr;
      importerNode["fileName"] = fileName.str();
      
#else
      auto impiGeometryNode = std::make_shared<ImpiSGNode>();
      impiGeometryNode->setName("impi_geometry");
      impiGeometryNode->setType("impi");

      // float values[8] = { 0,0,0,0,0,0,0,1 };
      // auto voxelArrayNode =
      //   std::make_shared<sg::DataArray1f>((float*)values,8,false);
      // voxelArrayNode->setName("voxel");
      // voxelArrayNode->setType("DataArray1f");
      // impiGeometryNode->add(voxelArrayNode);
      // impiGeometryNode->dims = 
#endif

      auto &impiMaterial = (*(*impiGeometryNode)["materialList"].nodeAs<sg::MaterialList>())[0];
      // auto &impiMaterial = (*impiGeometryNode)["material"];
      impiMaterial["Kd"] = vec3f(0.5f);
      impiMaterial["Ks"] = vec3f(0.1f);
      impiMaterial["Ns"] = 10.f;

      auto &lights = renderer["lights"];
      {
        auto &sun = lights.createChild("sun", "DirectionalLight");
        sun["color"] = vec3f(1.f,232.f/255.f,166.f/255.f);
        sun["direction"] = vec3f(0.462f,-1.f,-.1f);
        sun["intensity"] = 1.5f;

        auto &bounce = lights.createChild("bounce", "DirectionalLight");
        bounce["color"] = vec3f(127.f/255.f,178.f/255.f,255.f/255.f);
        bounce["direction"] = vec3f(-.93,-.54f,-.605f);
        bounce["intensity"] = 0.25f;
        
        auto &ambient = lights.createChild("ambient", "AmbientLight");
        ambient["intensity"] = 0.9f;
        ambient["color"] = vec3f(174.f/255.f,218.f/255.f,255.f/255.f);
      }
      
      
      // patchesInstance["model"].
      world.add(impiGeometryNode);

      ospray::ImGuiViewer window(renderer_ptr);

      auto &viewPort = window.viewPort;
      // XXX SG is too restrictive: OSPRay cameras accept non-normalized directions
      auto dir = normalize(viewPort.at - viewPort.from);
      renderer["camera"]["dir"] = dir;
      renderer["camera"]["pos"] = viewPort.from;
      renderer["camera"]["up"]  = viewPort.up;
      renderer["camera"]["fovy"] = viewPort.openingAngle;
      renderer["camera"]["apertureRadius"] = viewPort.apertureRadius;
      if (renderer["camera"].hasChild("focusdistance"))
        renderer["camera"]["focusdistance"] = length(viewPort.at - viewPort.from);

      window.create("OSPRay Example Viewer (module) App");

      ospray::imgui3D::run();
      return 0;
    }

  } // ::ospray::impi
} // ::ospray
