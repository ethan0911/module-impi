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

// static std::vector<vec3f> colors = {
//   vec3f(0.0, 0.000, 0.563),
//   vec3f(0.0, 0.000, 1.000),
//   vec3f(0.0, 1.000, 1.000),
//   vec3f(0.5, 1.000, 0.500),
//   vec3f(1.0, 1.000, 0.000),
//   vec3f(1.0, 0.000, 0.000),
//   vec3f(0.5, 0.000, 0.000),
// };
// static std::vector<float> opacities = { 1.f, 1.f };

static std::vector<float> colors = {
    0, 0, 1,
    0, 0.0130719, 1,
    0, 0.0261438, 1,
    0, 0.0392157, 1,
    0, 0.0522876, 1,
    0, 0.0653595, 1,
    0, 0.0784314, 1,
    0, 0.0915033, 1,
    0, 0.104575, 1,
    0, 0.117647, 1,
    0, 0.130719, 1,
    0, 0.143791, 1,
    0, 0.156863, 1,
    0, 0.169935, 1,
    0, 0.183007, 1,
    0, 0.196078, 1,
    0, 0.20915, 1,
    0, 0.222222, 1,
    0, 0.235294, 1,
    0, 0.248366, 1,
    0, 0.261438, 1,
    0, 0.27451, 1,
    0, 0.287582, 1,
    0, 0.300654, 1,
    0, 0.313726, 1,
    0, 0.326797, 1,
    0, 0.339869, 1,
    0, 0.352941, 1,
    0, 0.366013, 1,
    0, 0.379085, 1,
    0, 0.392157, 1,
    0, 0.405229, 1,
    0, 0.418301, 1,
    0, 0.431373, 1,
    0, 0.444444, 1,
    0, 0.457516, 1,
    0, 0.470588, 1,
    0, 0.48366, 1,
    0, 0.496732, 1,
    0, 0.509804, 1,
    0, 0.522876, 1,
    0, 0.535948, 1,
    0, 0.54902, 1,
    0, 0.562091, 1,
    0, 0.575163, 1,
    0, 0.588235, 1,
    0, 0.601307, 1,
    0, 0.614379, 1,
    0, 0.627451, 1,
    0, 0.640523, 1,
    0, 0.653595, 1,
    0, 0.666667, 1,
    0, 0.679739, 1,
    0, 0.69281, 1,
    0, 0.705882, 1,
    0, 0.718954, 1,
    0, 0.732026, 1,
    0, 0.745098, 1,
    0, 0.75817, 1,
    0, 0.771242, 1,
    0, 0.784314, 1,
    0, 0.797386, 1,
    0, 0.810458, 1,
    0, 0.823529, 1,
    0, 0.836601, 1,
    0, 0.849673, 1,
    0, 0.862745, 1,
    0, 0.875817, 1,
    0, 0.888889, 1,
    0, 0.901961, 1,
    0, 0.915033, 1,
    0, 0.928105, 1,
    0, 0.941176, 1,
    0, 0.954248, 1,
    0, 0.96732, 1,
    0, 0.980392, 1,
    0, 0.993464, 1,
    0.00653595, 1, 0.993464,
    0.0196078, 1, 0.980392,
    0.0326797, 1, 0.96732,
    0.0457516, 1, 0.954248,
    0.0588235, 1, 0.941176,
    0.0718954, 1, 0.928105,
    0.0849673, 1, 0.915033,
    0.0980392, 1, 0.901961,
    0.111111, 1, 0.888889,
    0.124183, 1, 0.875817,
    0.137255, 1, 0.862745,
    0.150327, 1, 0.849673,
    0.163399, 1, 0.836601,
    0.176471, 1, 0.823529,
    0.189542, 1, 0.810458,
    0.202614, 1, 0.797386,
    0.215686, 1, 0.784314,
    0.228758, 1, 0.771242,
    0.24183, 1, 0.75817,
    0.254902, 1, 0.745098,
    0.267974, 1, 0.732026,
    0.281046, 1, 0.718954,
    0.294118, 1, 0.705882,
    0.30719, 1, 0.69281,
    0.320262, 1, 0.679739,
    0.333333, 1, 0.666667,
    0.346405, 1, 0.653595,
    0.359477, 1, 0.640523,
    0.372549, 1, 0.627451,
    0.385621, 1, 0.614379,
    0.398693, 1, 0.601307,
    0.411765, 1, 0.588235,
    0.424837, 1, 0.575163,
    0.437909, 1, 0.562091,
    0.45098, 1, 0.54902,
    0.464052, 1, 0.535948,
    0.477124, 1, 0.522876,
    0.490196, 1, 0.509804,
    0.503268, 1, 0.496732,
    0.51634, 1, 0.48366,
    0.529412, 1, 0.470588,
    0.542484, 1, 0.457516,
    0.555556, 1, 0.444444,
    0.568627, 1, 0.431373,
    0.581699, 1, 0.418301,
    0.594771, 1, 0.405229,
    0.607843, 1, 0.392157,
    0.620915, 1, 0.379085,
    0.633987, 1, 0.366013,
    0.647059, 1, 0.352941,
    0.660131, 1, 0.339869,
    0.673203, 1, 0.326797,
    0.686275, 1, 0.313725,
    0.699346, 1, 0.300654,
    0.712418, 1, 0.287582,
    0.72549, 1, 0.27451,
    0.738562, 1, 0.261438,
    0.751634, 1, 0.248366,
    0.764706, 1, 0.235294,
    0.777778, 1, 0.222222,
    0.79085, 1, 0.20915,
    0.803922, 1, 0.196078,
    0.816993, 1, 0.183007,
    0.830065, 1, 0.169935,
    0.843137, 1, 0.156863,
    0.856209, 1, 0.143791,
    0.869281, 1, 0.130719,
    0.882353, 1, 0.117647,
    0.895425, 1, 0.104575,
    0.908497, 1, 0.0915033,
    0.921569, 1, 0.0784314,
    0.934641, 1, 0.0653595,
    0.947712, 1, 0.0522876,
    0.960784, 1, 0.0392157,
    0.973856, 1, 0.0261438,
    0.986928, 1, 0.0130719,
    1, 1, 0,
    1, 0.990196, 0,
    1, 0.980392, 0,
    1, 0.970588, 0,
    1, 0.960784, 0,
    1, 0.95098, 0,
    1, 0.941176, 0,
    1, 0.931373, 0,
    1, 0.921569, 0,
    1, 0.911765, 0,
    1, 0.901961, 0,
    1, 0.892157, 0,
    1, 0.882353, 0,
    1, 0.872549, 0,
    1, 0.862745, 0,
    1, 0.852941, 0,
    1, 0.843137, 0,
    1, 0.833333, 0,
    1, 0.823529, 0,
    1, 0.813725, 0,
    1, 0.803922, 0,
    1, 0.794118, 0,
    1, 0.784314, 0,
    1, 0.77451, 0,
    1, 0.764706, 0,
    1, 0.754902, 0,
    1, 0.745098, 0,
    1, 0.735294, 0,
    1, 0.72549, 0,
    1, 0.715686, 0,
    1, 0.705882, 0,
    1, 0.696078, 0,
    1, 0.686275, 0,
    1, 0.676471, 0,
    1, 0.666667, 0,
    1, 0.656863, 0,
    1, 0.647059, 0,
    1, 0.637255, 0,
    1, 0.627451, 0,
    1, 0.617647, 0,
    1, 0.607843, 0,
    1, 0.598039, 0,
    1, 0.588235, 0,
    1, 0.578431, 0,
    1, 0.568627, 0,
    1, 0.558823, 0,
    1, 0.549019, 0,
    1, 0.539216, 0,
    1, 0.529412, 0,
    1, 0.519608, 0,
    1, 0.509804, 0,
    1, 0.5, 0,
    1, 0.490196, 0,
    1, 0.480392, 0,
    1, 0.470588, 0,
    1, 0.460784, 0,
    1, 0.45098, 0,
    1, 0.441176, 0,
    1, 0.431372, 0,
    1, 0.421568, 0,
    1, 0.411765, 0,
    1, 0.401961, 0,
    1, 0.392157, 0,
    1, 0.382353, 0,
    1, 0.372549, 0,
    1, 0.362745, 0,
    1, 0.352941, 0,
    1, 0.343137, 0,
    1, 0.333333, 0,
    1, 0.323529, 0,
    1, 0.313725, 0,
    1, 0.303921, 0,
    1, 0.294118, 0,
    1, 0.284314, 0,
    1, 0.27451, 0,
    1, 0.264706, 0,
    1, 0.254902, 0,
    1, 0.245098, 0,
    1, 0.235294, 0,
    1, 0.22549, 0,
    1, 0.215686, 0,
    1, 0.205882, 0,
    1, 0.196078, 0,
    1, 0.186274, 0,
    1, 0.17647, 0,
    1, 0.166667, 0,
    1, 0.156863, 0,
    1, 0.147059, 0,
    1, 0.137255, 0,
    1, 0.127451, 0,
    1, 0.117647, 0,
    1, 0.107843, 0,
    1, 0.0980391, 0,
    1, 0.0882351, 0,
    1, 0.0784312, 0,
    1, 0.0686273, 0,
    1, 0.0588234, 0,
    1, 0.0490195, 0,
    1, 0.0392156, 0,
    1, 0.0294116, 0,
    1, 0.0196077, 0,
    1, 0.00980377, 0,
    1, 0, 0,
};
static std::vector<float> opacities = {
    0, 
    0.00392157, 
    0.00784314, 
    0.0117647, 
    0.0156863, 
    0.0196078, 
    0.0235294, 
    0.027451, 
    0.0313726, 
    0.0352941, 
    0.0392157, 
    0.0431373, 
    0.0470588, 
    0.0509804, 
    0.054902, 
    0.0588235, 
    0.0627451, 
    0.0666667, 
    0.0705882, 
    0.0745098, 
    0.0784314, 
    0.0823529, 
    0.0862745, 
    0.0901961, 
    0.0941177, 
    0.0980392, 
    0.101961, 
    0.105882, 
    0.109804, 
    0.113725, 
    0.117647, 
    0.121569, 
    0.12549, 
    0.129412, 
    0.133333, 
    0.137255, 
    0.141176, 
    0.145098, 
    0.14902, 
    0.152941, 
    0.156863, 
    0.160784, 
    0.164706, 
    0.168627, 
    0.172549, 
    0.176471, 
    0.180392, 
    0.184314, 
    0.188235, 
    0.192157, 
    0.196078, 
    0.2, 
    0.203922, 
    0.207843, 
    0.211765, 
    0.215686, 
    0.219608, 
    0.223529, 
    0.227451, 
    0.231373, 
    0.235294, 
    0.239216, 
    0.243137, 
    0.247059, 
    0.25098, 
    0.254902, 
    0.258824, 
    0.262745, 
    0.266667, 
    0.270588, 
    0.27451, 
    0.278431, 
    0.282353, 
    0.286275, 
    0.290196, 
    0.294118, 
    0.298039, 
    0.301961, 
    0.305882, 
    0.309804, 
    0.313726, 
    0.317647, 
    0.321569, 
    0.32549, 
    0.329412, 
    0.333333, 
    0.337255, 
    0.341176, 
    0.345098, 
    0.34902, 
    0.352941, 
    0.356863, 
    0.360784, 
    0.364706, 
    0.368627, 
    0.372549, 
    0.376471, 
    0.380392, 
    0.384314, 
    0.388235, 
    0.392157, 
    0.396078, 
    0.4, 
    0.403922, 
    0.407843, 
    0.411765, 
    0.415686, 
    0.419608, 
    0.423529, 
    0.427451, 
    0.431373, 
    0.435294, 
    0.439216, 
    0.443137, 
    0.447059, 
    0.45098, 
    0.454902, 
    0.458824, 
    0.462745, 
    0.466667, 
    0.470588, 
    0.47451, 
    0.478431, 
    0.482353, 
    0.486275, 
    0.490196, 
    0.494118, 
    0.498039, 
    0.496937, 
    0.490812, 
    0.484687, 
    0.478562, 
    0.472437, 
    0.466312, 
    0.460187, 
    0.454062, 
    0.447936, 
    0.441811, 
    0.435686, 
    0.429561, 
    0.423436, 
    0.417311, 
    0.411186, 
    0.405061, 
    0.398936, 
    0.39281, 
    0.386685, 
    0.38056, 
    0.374435, 
    0.36831, 
    0.362185, 
    0.35606, 
    0.349935, 
    0.343809, 
    0.337684, 
    0.331559, 
    0.325434, 
    0.319309, 
    0.313184, 
    0.307059, 
    0.300934, 
    0.294809, 
    0.288683, 
    0.282558, 
    0.276433, 
    0.270308, 
    0.264183, 
    0.258058, 
    0.251933, 
    0.245808, 
    0.239682, 
    0.233557, 
    0.227432, 
    0.221307, 
    0.215182, 
    0.209057, 
    0.202932, 
    0.196807, 
    0.190682, 
    0.184556, 
    0.178431, 
    0.172306, 
    0.166181, 
    0.160056, 
    0.153931, 
    0.147806, 
    0.141681, 
    0.135556, 
    0.12943, 
    0.123305, 
    0.11718, 
    0.111055, 
    0.10493, 
    0.0988047, 
    0.0926796, 
    0.0865545, 
    0.0804294, 
    0.0743043, 
    0.0681792, 
    0.062054, 
    0.0559289, 
    0.0498038, 
    0.0436787, 
    0.0375536, 
    0.0314285, 
    0.0253033, 
    0.0191782, 
    0.0130531, 
    0.006928, 
    0.000802875, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
};

int main(int ac, const char** av)
{
  //-----------------------------------------------------
  // Program Initialization
  //----------------------------------------------------- 
#ifdef __unix__
  // check hostname
  char hname[200];
  gethostname(hname, 200);
  std::cout << "#osp: on host >> " << hname << " <<" << std::endl;;
#endif
  int init_error = ospInit(&ac, av);
  //-----------------------------------------------------
  // Master Rank Code (worker nodes will not reach here)
  //-----------------------------------------------------

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

#if USE_VIEWER
  int window = ospray::viewer::Init(ac, av, imgSize.x, imgSize.y);
#endif

  //-----------------------------------------------------
  // Create ospray context
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
  // Create ospray objects
  //-----------------------------------------------------  

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
  OSPData colorsData = ospNewData(colors.size() / 3, OSP_FLOAT3,
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

  ospray::viewer::Handler(camera, (const osp::vec3f&)vp, (const osp::vec3f&)vu, (const osp::vec3f&)vi);
  ospray::viewer::Handler(transferFcn, amrVolume->Range().x, amrVolume->Range().y);
  ospray::viewer::Handler(world, renderer);
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
