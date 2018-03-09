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
#include "ospray/volume/amr/AMRVolume.h"
#include "ospcommon/utility/getEnvVar.h"
#include "exampleViewer/widgets/imguiViewer.h"
#include "impiCommandLine.h"

/*! _everything_ in the ospray core universe should _always_ be in the
  'ospray' namespace. */
namespace ospray {

  /*! though not required, it is good practice to put any module into
    its own namespace (isnide of ospray:: ). Unlike for the naming of
    library and init function, the naming for this namespace doesn't
    particularlly matter. E.g., 'impi', 'module_blp',
    'bilinar_patch' etc would all work equally well. */
  namespace impi {

    struct clTransform
    {
      vec3f translate{0, 0, 0};
      vec3f scale{.5f, .5f, .5f};
      vec3f rotation{0, 0, 0};
    };

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
                      const unsigned char *binBasePtr) override
      {
      }

      void postCommit(sg::RenderContext &ctx)
      {
        auto materialListNode    = child("materialList").nodeAs<sg::MaterialList>();
        const auto &materialList = materialListNode->nodes;
        if (!materialList.empty()) {
          std::vector<OSPObject> mats;
          for (auto mat : materialList) {
            auto m = mat->valueAs<OSPObject>();
            if (m)
              mats.push_back(m);
          }
          auto ospMaterialList =
              ospNewData(mats.size(), OSP_OBJECT, mats.data());
          ospCommit(ospMaterialList);
          ospSetData(valueAs<OSPObject>(), "materialList", ospMaterialList);
        }

        Geometry::postCommit(ctx);
      }
    };

    // use ospcommon for vec3f etc
    using namespace ospcommon;

    std::shared_ptr<sg::Importer> importObject2World(sg::Node& renderer, std::string fileName)
    {
      auto& world = renderer["world"];
      auto importerNode_ptr =
          sg::createNode(fileName, "Importer")->nodeAs<sg::Importer>();

      auto &importerNode = *importerNode_ptr;
      importerNode["fileName"] = fileName;

      clTransform cltransform;
      auto &transform =
          world.createChild("transform_" + fileName, "Transform");
      transform["scale"]    = cltransform.scale;
      transform["rotation"] = cltransform.rotation;

      transform.add(importerNode_ptr);
      renderer.traverse("verify");
      renderer.traverse("commit");
      auto bounds   = importerNode_ptr->computeBounds();
      auto size     = bounds.upper - bounds.lower;
      float maxSize = max(max(size.x, size.y), size.z);
      if (!std::isfinite(maxSize))
        maxSize = 0.f;  // FIXME: why is maxSize = NaN in some cases?!
      vec3f offset          = {maxSize * 1.3f, maxSize * 1.3f, maxSize * 1.3f};
      transform["position"] = cltransform.translate + offset;
      
      return importerNode_ptr;
    }

    //----------------------------------------------------------------------------
    // Here we add all the dataset specific ha
    //----------------------------------------------------------------------------
    struct DataState {
      float isoValue = 0.f;
      vec3f sunDir = vec3f(-1.f,0.679f,-0.754f);
      vec3f disDir = vec3f(.372f,.416f,-0.605f);
    } dataState;
    std::shared_ptr<sg::Node> DataParser(const std::string& dataString,
					 const CommandLine &args, 
					 sg::Node &renderer, sg::Node &world, 
					 std::shared_ptr<sg::Node> impiGeometryNode)
    {
      if (dataString == "landingGear") {
	// Check Arguments
	if (args.inputFiles.size() != 2) {
	  throw std::runtime_error("missing input arguments for 'landingGear' dataset");
	}
	// Setup Data States
	dataState.sunDir = vec3f(.783f, -1.f, -0.086f);
	dataState.disDir = vec3f(.337f, .416f, -0.605f);
	dataState.isoValue = 99000.0f;
	// Set Geometry (qi: what it is doing here?) 
	auto model = sg::createNode("Impl_model", "Model");
	model->add(impiGeometryNode);
	auto objInstance = sg::createNode("instance", "Instance");
	objInstance->setChild("model", model);
	model->setParent(objInstance);
	auto landingGearImportNodePtr = importObject2World(renderer, args.inputFiles[1]);
	auto& instance = (*landingGearImportNodePtr)["instance"];
	instance.child("position").setValue(vec3f(-61.61,-61.6,-93.4));
	instance.child("scale").setValue(vec3f(2.f,2.f,2.f));	
	return objInstance;
      }
      else if (dataString == "chombo") {
	dataState.isoValue = 0.7f;
	return impiGeometryNode;
      }
      else if (dataString == "chombo") {
	dataState.isoValue = 0.0f;
	return impiGeometryNode;
      }
      else {
	return impiGeometryNode;
      }
    }
    //----------------------------------------------------------------------------

    extern "C" int main(int ac, const char **av)
    {
      //-----------------------------------------------------
      // Program Initialization
      //-----------------------------------------------------
      int init_error = ospInit(&ac, av);
      
      //-----------------------------------------------------
      // Master Rank Code (worker nodes will not reach here)
      //-----------------------------------------------------
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
      loadLibrary("ospray_sg"); // access/load symbols/sg::Nodes dynamically
      ospLoadModule("impi");
      ospray::imgui3D::init(&ac,av);

      //-----------------------------------------------------
      // parse the commandline;
      // complain about anything we do not recognize
      //-----------------------------------------------------
      std::string dataString = 
	ospcommon::utility::getEnvVar<std::string>("IMPI_AMR_DATA").value_or("cosmos");     
      CommandLine args(ac, av);
      
      //-----------------------------------------------------
      // Setup Scene Graph
      //-----------------------------------------------------
      auto renderer_ptr = sg::createNode("renderer", "Renderer");
      auto &renderer = *renderer_ptr;

      auto &win_size = ospray::imgui3D::ImGui3DWidget::defaultInitSize;
      renderer["frameBuffer"]["size"] = win_size;

      renderer["rendererType"] = std::string("scivis");
      auto &world              = renderer["world"];

      if (!args.inputFiles.empty()) {
	
        auto importerNodePtr = importObject2World(renderer, args.inputFiles[0]);

        auto amrVolSGNodePtr = importerNodePtr->childByType("AMRVolume"); // foundNode->second;
        (*amrVolSGNodePtr)["visible"] = false; // AMR volume transfer function is buggy ?

        auto impiGeometryNode = std::make_shared<ImpiSGNode>();
        impiGeometryNode->setName("impi_geometry");
        impiGeometryNode->setType("impi");

	auto objNode = DataParser(dataString, args, renderer, world, impiGeometryNode);

        impiGeometryNode->createChild("isoValue", "float", dataState.isoValue);
        impiGeometryNode->createChild("amrDataPtr", "OSPObject", 
				      (OSPObject)amrVolSGNodePtr->valueAs<OSPVolume>());

        auto &impiMaterial = (*(*impiGeometryNode)["materialList"].nodeAs<sg::MaterialList>())[0];
        impiMaterial["Kd"] = vec3f(0.5f);
        impiMaterial["Ks"] = vec3f(0.1f);
        impiMaterial["Ns"] = 10.f;
	
	world.add(objNode);

      }

      //-----------------------------------------------------
      // Setup Lights
      //-----------------------------------------------------
      auto &lights = renderer["lights"];
      {
        auto &sun = lights.createChild("sun", "DirectionalLight");
        sun["color"] = vec3f(1.f,255.f/255.f,255.f/255.f);
        sun["direction"] = dataState.sunDir;
        sun["intensity"] = 1.5f;

        auto &bounce = lights.createChild("bounce", "DirectionalLight");
        bounce["color"] = vec3f(127.f/255.f,178.f/255.f,255.f/255.f);
        bounce["direction"] = dataState.disDir;
        bounce["intensity"] = 0.25f;

        auto &ambient = lights.createChild("ambient", "AmbientLight");
        ambient["intensity"] = 0.9f;
        ambient["color"] = vec3f(174.f/255.f,218.f/255.f,255.f/255.f);
      }

      //-----------------------------------------------------
      // Render with OpenGL
      //-----------------------------------------------------
      ospray::ImGuiViewer window(renderer_ptr);
      auto &viewPort = window.viewPort;
      // XXX SG is too restrictive: OSPRay cameras accept non-normalized directions
      auto dir = normalize(viewPort.at - viewPort.from);
      renderer["camera"]["dir"] = dir;
      renderer["camera"]["pos"] = viewPort.from;
      renderer["camera"]["up"]  = viewPort.up;
      renderer["camera"]["fovy"] = viewPort.openingAngle;
      renderer["camera"]["apertureRadius"] = viewPort.apertureRadius;
      if (renderer["camera"].hasChild("focusdistance")) {
        renderer["camera"]["focusdistance"] = length(viewPort.at - viewPort.from);
      }
      window.create("OSPRay Example Viewer (module) App");
      ospray::imgui3D::run();

      return 0;
    }

  } // ::ospray::impi
} // ::ospray
