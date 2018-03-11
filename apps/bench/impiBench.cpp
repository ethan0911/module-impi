// ======================================================================== //
// Copyright SCI Institute, University of Utah, 2018
// ======================================================================== //
// 
// Useful Parameters
// landingGear: 
//    -vp 16.286070 16.446814 0.245150
//    -vu -0.000000 -0.000000 -1.000000 
//    -vi 16.430407 16.157639 0.353916 
//    -scale 1 1 1 -translate 15.995 16 0.1
//
// ======================================================================== //

#include "ospray/ospray.h"

#include "ospcommon/vec.h"
#include "ospcommon/box.h"
#include "ospcommon/range.h"
#include "ospcommon/LinearSpace.h"
#include "ospcommon/AffineSpace.h"

#include "impiHelper.h"
#include "impiReader.h"
#include "loader/meshloader.h"

#ifdef __unix__
# include <unistd.h>
#endif

#if USE_VIEWER
#include "opengl/viewer.h"
#endif

using namespace ospcommon;

static bool showVolume{false};
static bool showObject{false};
static enum {IMPI, NORMAL} isoMode;
static std::string rendererName = "scivis";
static std::string outputImageName = "result";

static std::vector<float> isoValues(1, 0.0f);
static vec3f isoScale{1.f, 1.f, 1.f};
static vec3f isoTranslate{.0f,.0f,.0f};

static vec3f vp{43.2,44.9,-57.6};
static vec3f vu{0,1,0};
static vec3f vi{0,0,0};
static vec3f sunDir{-1.f,0.679f,-0.754f};
static vec3f disDir{.372f,.416f,-0.605f};
static vec2i imgSize{1024, 768};
static vec2i numFrames{1/* skipped */, 20/* measure */};
static affine3f Identity(vec3f(1,0,0), vec3f(0,1,0), vec3f(0,0,1), vec3f(0,0,0));
static std::vector<vec3f> colors = {
  vec3f(0.0, 0.000, 0.563),
  vec3f(0.0, 0.000, 1.000),
  vec3f(0.0, 1.000, 1.000),
  vec3f(0.5, 1.000, 0.500),
  vec3f(1.0, 1.000, 0.000),
  vec3f(1.0, 0.000, 0.000),
  vec3f(0.5, 0.000, 0.000),
};
static std::vector<float> opacities = { 1.f, 1.f };

int main(int ac, const char** av)
{
  //-----------------------------------------------------
  // Program Initialization
  //----------------------------------------------------- 
  // check hostname
#ifdef __unix__
  char hname[200];
  gethostname(hname, 200);
  std::cout << "#osp: on host >> " << hname << " <<" << std::endl;;
#endif
  int init_error = ospInit(&ac, av);
#if USE_VIEWER
  int window = ospray::viewer::Init(ac, av);
#endif

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
  std::vector<std::string> inputFiles;
  std::string              inputMesh;
  for (int i = 1; i < ac; ++i) {
    std::string str(av[i]);
    if (str == "-o") {
      outputImageName = av[++i];
    }
    else if (str == "-renderer") {
      rendererName = av[++i];
    }
    else if (str == "-iso" || str == "-isoValue") {
      ospray::impi::Parse<1>(ac, av, i, isoValues[0]);
    }
    else if (str == "-isos" || str == "-isoValues") {
      try {
	int n = 0;
	ospray::impi::Parse<1>(ac, av, i, n);
	isoValues.resize(n);
	for (int j = 0; j < n; ++j) {
	  float v = 0.0f;
	  ospray::impi::Parse<1>(ac, av, i, v);
	  isoValues[j] = v;
	}
      } catch (const std::runtime_error& e) {
	throw std::runtime_error(std::string(e.what())+
				 " usage: -isos "
				 "<# of values> "
				 "<value list>");
      }
    }
    else if (str == "-translate") {
      ospray::impi::Parse<3>(ac, av, i, isoTranslate);
    }
    else if (str == "-scale") {
      ospray::impi::Parse<3>(ac, av, i, isoScale);
    }
    else if (str == "-fb") {
      ospray::impi::Parse<2>(ac, av, i, imgSize);
    }
    else if (str == "-vp") { 
      ospray::impi::Parse<3>(ac, av, i, vp);
    }
    else if (str == "-vi") {
      ospray::impi::Parse<3>(ac, av, i, vi);
    }
    else if (str == "-vu") { 
      ospray::impi::Parse<3>(ac, av, i, vu);
    }    
    else if (str == "-sun") { 
      ospray::impi::Parse<3>(ac, av, i, sunDir);
    }    
    else if (str == "-dis") { 
      ospray::impi::Parse<3>(ac, av, i, disDir);
    }    
    else if (str == "-volume") { 
      showVolume = true;
    }
    else if (str == "-object") { 
      showObject = true;
      inputMesh = av[++i];
    }
    else if (str == "-use-builtin-isosurface") { 
      isoMode = NORMAL;
    }   
    else if (str == "-frames") {
      try {
	ospray::impi::Parse<2>(ac, av, i, numFrames);
      } catch (const std::runtime_error& e) {
	throw std::runtime_error(std::string(e.what())+
				 " usage: -frames "
				 "<# of warmup frames> "
				 "<# of benchmark frames>");
      }
    }
    else if (str[0] == '-') {
      throw std::runtime_error("unknown argument: " + str);
    }
    else {
      inputFiles.push_back(av[i]);
    }
  }
  if (inputFiles.empty()) { throw std::runtime_error("missing input file"); }
  if (inputFiles.size() > 1) {
    for (auto& s : inputFiles) std::cout << s << std::endl;
    throw std::runtime_error("too many input file"); 
  }
  
  // create world and renderer
  OSPModel world = ospNewModel();
  OSPRenderer renderer = ospNewRenderer(rendererName.c_str());
  std::cout << "#osp:bench: OSPRay renderer: " << rendererName << std::endl; 
  if (!renderer) {
    throw std::runtime_error("invalid renderer name: " + rendererName);
  }

  // load amr volume
  auto amrVolume = ospray::ParseOSP::loadOSP(inputFiles[0]);

  // setup trasnfer function
  OSPData colorsData = ospNewData(colors.size(), OSP_FLOAT3,
				  colors.data());
  ospCommit(colorsData);
  OSPData opacitiesData = ospNewData(opacities.size(), OSP_FLOAT, 
				     opacities.data());
  ospCommit(opacitiesData);
  OSPTransferFunction transferFcn = ospNewTransferFunction("piecewise_linear");
  ospSetData(transferFcn, "colors",    colorsData);
  ospSetData(transferFcn, "opacities", opacitiesData);
  ospSetVec2f(transferFcn, "valueRange", 
  	      (const osp::vec2f&)amrVolume->Range());
  ospCommit(transferFcn);
  ospRelease(colorsData);
  ospRelease(opacitiesData);

  // setup volume
  OSPVolume volume = amrVolume->Create(transferFcn);
  if (showVolume) {
    ospAddVolume(world, volume);
  }

  // setup isosurfaces
  OSPModel local = ospNewModel();
  // Node: all isosurfaces will share one material for now
  OSPMaterial mtl = ospNewMaterial(renderer, "OBJMaterial");
  ospSetVec3f(mtl, "Kd", osp::vec3f{0.5f, 0.5f, 0.5f});
  ospSetVec3f(mtl, "Ks", osp::vec3f{0.1f, 0.1f, 0.1f});
  ospSet1f(mtl, "Ns", 10.f);
  ospCommit(mtl);
  switch (isoMode) {
  case NORMAL:
    // --> normal isosurface
    {
      OSPGeometry niso = ospNewGeometry("isosurfaces");
      OSPData niso_values = ospNewData(isoValues.size(), 
				       OSP_FLOAT, 
				       isoValues.data());
      ospSetData(niso, "isovalues", niso_values);
      ospSetObject(niso, "volume", volume);
      ospSetMaterial(niso, mtl);
      ospCommit(niso);
      ospAddGeometry(local, niso);
      ospCommit(local);    
    }
    break;
  case IMPI:
    // --> implicit isosurface
    {
      // Note: because there is no naive multi-iso surface support,
      //       we build multiple iso-geometries here
      for (auto& v : isoValues) {
	OSPGeometry iiso = ospNewGeometry("impi"); 
	ospSet1f(iiso, "isoValue", v);
	ospSetObject(iiso, "amrDataPtr", volume);
	ospSetMaterial(iiso, mtl);
	ospCommit(iiso);
	ospAddGeometry(local, iiso);
      }
    }
    break;
  default:
    throw std::runtime_error("wrong ISO-Mode, this shouldn't happen");
  }
  ospCommit(local);  
  OSPGeometry isoinstance = 
    ospNewInstance(local, (const osp::affine3f &)Identity);
  ospAddGeometry(world, isoinstance);

  // setup object
  Mesh mesh;
  affine3f transform = 
    affine3f::translate(isoTranslate) * 
    affine3f::scale(isoScale);
  if (showObject) {
    mesh.LoadFromFileObj(inputMesh.c_str(), false);
    mesh.SetTransform(transform);
    mesh.AddToModel(world, renderer, mtl);
  }

  // setup camera
  OSPCamera camera = ospNewCamera("perspective");
  vec3f vd = vi - vp;
  ospSetVec3f(camera, "pos", (const osp::vec3f&)vp);
  ospSetVec3f(camera, "dir", (const osp::vec3f&)vd);
  ospSetVec3f(camera, "up",  (const osp::vec3f&)vu);
  ospSet1f(camera, "aspect", imgSize.x / (float)imgSize.y);
  ospSet1f(camera, "fovy", 60.f);
  ospCommit(camera);

  // setup lighting
  OSPLight d_light = ospNewLight(renderer, "DirectionalLight");
  ospSet1f(d_light, "intensity", 0.25f);
  ospSet1f(d_light, "angularDiameter", 0.53f);
  ospSetVec3f(d_light, "color", 
	      osp::vec3f{127.f/255.f,178.f/255.f,255.f/255.f});
  ospSetVec3f(d_light, "direction", (const osp::vec3f&)disDir);
  ospCommit(d_light);
  OSPLight s_light = ospNewLight(renderer, "DirectionalLight");
  ospSet1f(s_light, "intensity", 1.50f);
  ospSet1f(s_light, "angularDiameter", 0.53f);
  ospSetVec3f(s_light, "color", osp::vec3f{1.f,1.f,1.f});  
  ospSetVec3f(s_light, "direction", (const osp::vec3f&)sunDir);
  ospCommit(s_light);
  OSPLight a_light = ospNewLight(renderer, "AmbientLight");
  ospSet1f(a_light, "intensity", 0.90f);
  ospSetVec3f(a_light, "color", 
	      osp::vec3f{174.f/255.f,218.f/255.f,255.f/255.f});
  ospCommit(a_light);
  std::vector<OSPLight> light_list { a_light, d_light, s_light };
  OSPData lights = ospNewData(light_list.size(), OSP_OBJECT, 
			      light_list.data());
  ospCommit(lights);

  // setup world & renderer
  ospCommit(world); 
  ospSetVec3f(renderer, "bgColor", 
	      osp::vec3f{230.f/255.f, 230.f/255.f, 230.f/255.f});
  ospSetData(renderer, "lights", lights);
  ospSetObject(renderer, "model", world);
  ospSetObject(renderer, "camera", camera);
  ospSet1i(renderer, "shadowEnabled", 1);
  ospSet1i(renderer, "oneSidedLighting", 1);
  ospSet1i(renderer, "maxDepth", 5);
  ospSet1i(renderer, "spp", 1);
  ospSet1i(renderer, "autoEpsilon", 1);
  ospSet1i(renderer, "aoSamples", 1);
  ospSet1i(renderer, "aoTransparencyEnabled", 1);
  ospSet1f(renderer, "aoDistance", 10000.0f);
  ospSet1f(renderer, "epsilon", 0.001f);
  ospSet1f(renderer, "minContribution", 0.001f);
  ospCommit(renderer);

#if USE_VIEWER

  ospray::viewer::Handler(camera);
  ospray::viewer::Handler(transferFcn, amrVolume->Range().x, amrVolume->Range().y);
  ospray::viewer::Handler(renderer);
  ospray::viewer::Render(window);

#else

  // setup framebuffer
  OSPFrameBuffer fb = ospNewFrameBuffer((const osp::vec2i&)imgSize, 
					OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
  ospFrameBufferClear(fb, OSP_FB_COLOR | OSP_FB_ACCUM);

  // render 10 more frames, which are accumulated to result in a better converged image
  std::cout << "#osp:bench: start warmups for " << numFrames.x << " frames" << std::endl;
  for (int frames = 0; frames < numFrames.x; frames++) {   // skip some frames to warmup
    ospRenderFrame(fb, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
  }
  std::cout << "#osp:bench: done warmups" << std::endl;
  std::cout << "#osp:bench: start benchmarking for " << numFrames.y << " frames" << std::endl;
  auto t = ospray::impi::Time();
  for (int frames = 0; frames < numFrames.y; frames++) {
    ospRenderFrame(fb, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
  }
  auto et = ospray::impi::Time(t);
  std::cout << "#osp:bench: done benchmarking" << std::endl;
  std::cout << "#osp:bench: average framerate: " << numFrames.y/et << std::endl; 

  // save frame
  const uint32_t * buffer = (uint32_t*)ospMapFrameBuffer(fb, OSP_FB_COLOR);
  ospray::impi::writePPM(outputImageName + ".ppm", imgSize.x, imgSize.y, buffer);
  ospUnmapFrameBuffer(buffer, fb);

#endif

  // done
  std::cout << "#osp:bench: done benchmarking" << std::endl;
  return 0;
}
