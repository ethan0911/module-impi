// ======================================================================== //
// Copyright SCI Institute, University of Utah, 2018
// ======================================================================== //

#include "ospray/ospray.h"
#include "impiHelper.h"
#include "impiReader.h"
#include "impiCommandLine.h"

using namespace ospcommon;

static float isoValue = 0.0f;
static vec2i imgSize{1024, 768};
static vec3f vp{43.2,44.9,-57.6};
static vec3f vd{-0.534,-0.456,0.712};;
static vec3f vu{0,1,0};
static std::vector<vec3f> colors = {
  vec3f(0.0, 0.000, 0.563),
  vec3f(0.0, 0.000, 1.000),
  vec3f(0.0, 1.000, 1.000),
  vec3f(0.5, 1.000, 0.500),
  vec3f(1.0, 1.000, 0.000),
  vec3f(1.0, 0.000, 0.000),
  vec3f(0.5, 0.000, 0.000),
};
static std::vector<float> opacities = { 0.01f, 0.01f };

int main(int ac, const char** av)
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
  ospLoadModule("impi");

  //-----------------------------------------------------
  // parse the commandline;
  // complain about anything we do not recognize
  //-----------------------------------------------------
  ospray::impi::CommandLine args;
  try {
    args.Parse(ac, av);
  } catch (const std::runtime_error& error) {
    for (int i = 1; i < ac; ++i) {
      std::string str(av[i]);
      if (str == "-iso" || str == "-isoValue") {
	ospray::impi::Parse<1>(ac, av, i, isoValue);
      }
      else if (str == "-fb") {
	ospray::impi::Parse<2>(ac, av, i, imgSize);
      }
      else if (str == "-vp") { 
	ospray::impi::Parse<3>(ac, av, i, vp);
      }
      else if (str == "-vd") { 
	ospray::impi::Parse<3>(ac, av, i, vd);
      }
      else if (str == "-vu") { 
	ospray::impi::Parse<3>(ac, av, i, vu);
      }    
    }
  }
  if (args.inputFiles.empty()) { throw std::runtime_error("missing input file"); }

  // create world and renderer
  OSPModel world = ospNewModel();
  OSPRenderer renderer = ospNewRenderer("scivis");

  // setup trasnfer function
  OSPData colorsData = ospNewData(colors.size(), OSP_FLOAT3, colors.data());
  ospCommit(colorsData);
  OSPData opacitiesData = ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
  ospCommit(opacitiesData);
  OSPTransferFunction transferFcn = ospNewTransferFunction("piecewise_linear");
  ospSetData(transferFcn, "colors",    colorsData);
  ospSetData(transferFcn, "opacities", opacitiesData);
  ospCommit(transferFcn);
  ospRelease(colorsData);
  ospRelease(opacitiesData);

  // setup volume
  auto amrVolume = ospray::ParseOSP::loadOSP(args.inputFiles[0]);
  OSPVolume volume = amrVolume->Create(transferFcn);
  ospAddVolume(world, volume);

  // setup isosurface
  OSPGeometry isosurface = ospNewGeometry("impi"); 
  ospSet1f(isosurface, "isoValue", isoValue);
  ospSetObject(isosurface, "amrDataPtr", volume);
  OSPMaterial mtl = ospNewMaterial(renderer, "OBJMaterial");
  ospSetVec3f(mtl, "Kd", osp::vec3f{0.5f, 0.5f, 0.5f});
  ospSetVec3f(mtl, "Ks", osp::vec3f{0.1f, 0.1f, 0.1f});
  ospSet1f(mtl, "Ns", 10.f);
  ospCommit(mtl);
  ospSetMaterial(isosurface, mtl);
  ospCommit(isosurface);
  ospAddGeometry(world, isosurface);

  // setup camera
  OSPCamera camera = ospNewCamera("perspective");
  ospSetVec3f(camera, "pos", (const osp::vec3f&)vp);
  ospSetVec3f(camera, "dir", (const osp::vec3f&)vd);
  ospSetVec3f(camera, "up",  (const osp::vec3f&)vu);
  ospSet1f(camera, "aspect", imgSize.x / (float)imgSize.y);
  ospSet1f(camera, "fovy", 60.f);
  ospCommit(camera);

  // setup lighting
  OSPLight a_light = ospNewLight(renderer, "AmbientLight");
  ospSet1f(a_light, "intensity", 0.90f);
  ospSetVec3f(a_light, "color", osp::vec3f{174.f/255.f,218.f/255.f,255.f/255.f});
  ospCommit(a_light);
  OSPLight d_light = ospNewLight(renderer, "DirectionalLight");
  ospSet1f(d_light, "intensity", 0.25f);
  ospSetVec3f(d_light, "color", osp::vec3f{127.f/255.f,178.f/255.f,255.f/255.f});
  ospSetVec3f(d_light, "direction", osp::vec3f{.372f,.416f,-0.605f});
  ospCommit(d_light);
  OSPLight s_light = ospNewLight(renderer, "DirectionalLight");
  ospSet1f(s_light, "intensity", 1.50f);
  ospSetVec3f(s_light, "color", osp::vec3f{1.f,255.f/255.f,255.f/255.f});
  ospSetVec3f(s_light, "direction", osp::vec3f{-1.f,0.679f,-0.754f});
  ospCommit(s_light);
  std::vector<OSPLight> light_list { a_light, d_light, s_light };
  OSPData lights = ospNewData(light_list.size(), OSP_OBJECT, light_list.data());
  ospCommit(lights);

  // setup world & renderer
  ospCommit(world); 
  ospSetVec3f(renderer, "bgColor", osp::vec3f{0.5f, 0.5f, 0.5f});
  ospSetData(renderer, "lights", lights);
  ospSetObject(renderer, "model", world);
  ospSetObject(renderer, "camera", camera);
  ospSet1i(renderer, "shadowEnabled", 0);
  ospSet1i(renderer, "oneSidedLighting", 0);
  ospCommit(renderer);

  // setup framebuffer
  OSPFrameBuffer fb = ospNewFrameBuffer((const osp::vec2i&)imgSize, 
					OSP_FB_SRGBA,
					OSP_FB_COLOR | 
					OSP_FB_ACCUM);
  ospFrameBufferClear(fb, OSP_FB_COLOR | OSP_FB_ACCUM);

  // render 10 more frames, which are accumulated to result in a better converged image
  for (int frames = 0; frames < 10; frames++) {
    ospRenderFrame(fb, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
  }

  // save frame
  const uint32_t * buffer = (uint32_t*)ospMapFrameBuffer(fb, OSP_FB_COLOR);
  writePPM("result.ppm", imgSize.x, imgSize.y, buffer);
  ospUnmapFrameBuffer(buffer, fb);

  return 0;
}
