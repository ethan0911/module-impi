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

struct ISO {
  float v = 0.0f;
  vec3f c = vec3f(0.5f, 0.5f, 0.5f);
  OSPMaterial mtl;
  OSPGeometry geo;
};
static std::vector<ISO> isoValues(1);

static vec3f objScale{1.f, 1.f, 1.f};
static vec3f objTranslate{.0f,.0f,.0f};

static vec3f vp{43.2,44.9,-57.6};
static vec3f vu{0,1,0};
static vec3f vi{0,0,0};
static vec3f sunDir{-1.f,0.679f,-0.754f};
static vec3f disDir{.372f,.416f,-0.605f};
static vec2i imgSize{1024, 768};
static vec2i numFrames{1/* skipped */, 20/* measure */};
static affine3f Identity(vec3f(1,0,0), vec3f(0,1,0), vec3f(0,0,1), vec3f(0,0,0));
static std::vector<float> colors = {
    0, 0, 0,
    0, 0.00755413, 0.0189916,
    0, 0.0151083, 0.0379831,
    0, 0.0226624, 0.0569747,
    0, 0.0302165, 0.0759662,
    0, 0.0377707, 0.0949578,
    0, 0.0453248, 0.113949,
    0, 0.0528789, 0.132941,
    0, 0.0604331, 0.151932,
    0, 0.0679872, 0.170924,
    0, 0.0755413, 0.189916,
    0, 0.0830955, 0.208907,
    0, 0.0906496, 0.227899,
    0, 0.0982037, 0.24689,
    0, 0.105758, 0.265882,
    0, 0.113312, 0.284873,
    0, 0.120771, 0.303548,
    0, 0.126807, 0.317471,
    0, 0.132843, 0.331394,
    0, 0.138878, 0.345317,
    0, 0.144914, 0.35924,
    0, 0.150949, 0.373163,
    0, 0.156985, 0.387086,
    0, 0.163021, 0.401009,
    0, 0.169056, 0.414932,
    0, 0.175092, 0.428855,
    0, 0.181128, 0.442778,
    0, 0.187163, 0.456701,
    0, 0.193199, 0.470624,
    0, 0.199235, 0.484547,
    0, 0.20527, 0.498469,
    0, 0.211306, 0.512392,
    0.000433363, 0.217594, 0.525633,
    0.00390021, 0.225653, 0.534099,
    0.00736706, 0.233712, 0.542564,
    0.0108339, 0.24177, 0.55103,
    0.0143008, 0.249829, 0.559495,
    0.0177676, 0.257888, 0.567961,
    0.0212345, 0.265946, 0.576427,
    0.0247013, 0.274005, 0.584892,
    0.0281681, 0.282064, 0.593358,
    0.031635, 0.290122, 0.601823,
    0.0351018, 0.298181, 0.610289,
    0.0385687, 0.30624, 0.618755,
    0.0420355, 0.314298, 0.62722,
    0.0455024, 0.322357, 0.635686,
    0.0489692, 0.330416, 0.644151,
    0.0524361, 0.338474, 0.652617,
    0.0561094, 0.346758, 0.66021,
    0.0606773, 0.356017, 0.664025,
    0.0652452, 0.365277, 0.667839,
    0.0698131, 0.374536, 0.671653,
    0.0743811, 0.383795, 0.675468,
    0.078949, 0.393055, 0.679282,
    0.0835169, 0.402314, 0.683097,
    0.0880848, 0.411573, 0.686911,
    0.0926527, 0.420833, 0.690725,
    0.0972206, 0.430092, 0.69454,
    0.101789, 0.439351, 0.698354,
    0.106356, 0.44861, 0.702169,
    0.110924, 0.45787, 0.705983,
    0.115492, 0.467129, 0.709797,
    0.12006, 0.476388, 0.713612,
    0.124628, 0.485648, 0.717426,
    0.129009, 0.494925, 0.721413,
    0.13283, 0.504256, 0.725919,
    0.136651, 0.513587, 0.730425,
    0.140472, 0.522918, 0.73493,
    0.144293, 0.532249, 0.739436,
    0.148115, 0.54158, 0.743942,
    0.151936, 0.550911, 0.748447,
    0.155757, 0.560242, 0.752953,
    0.159578, 0.569573, 0.757459,
    0.163399, 0.578904, 0.761964,
    0.16722, 0.588236, 0.76647,
    0.171041, 0.597567, 0.770976,
    0.174862, 0.606898, 0.775481,
    0.178683, 0.616229, 0.779987,
    0.182504, 0.62556, 0.784493,
    0.186325, 0.634891, 0.788998,
    0.191672, 0.644122, 0.793691,
    0.200376, 0.653135, 0.798794,
    0.20908, 0.662147, 0.803897,
    0.217784, 0.671159, 0.809,
    0.226488, 0.680171, 0.814103,
    0.235192, 0.689184, 0.819206,
    0.243896, 0.698196, 0.824309,
    0.2526, 0.707208, 0.829412,
    0.261304, 0.71622, 0.834515,
    0.270008, 0.725233, 0.839618,
    0.278712, 0.734245, 0.844721,
    0.287416, 0.743257, 0.849824,
    0.29612, 0.75227, 0.854927,
    0.304824, 0.761282, 0.86003,
    0.313528, 0.770294, 0.865134,
    0.322232, 0.779306, 0.870237,
    0.334274, 0.787462, 0.874888,
    0.351878, 0.79419, 0.878785,
    0.369482, 0.800918, 0.882683,
    0.387087, 0.807646, 0.886581,
    0.404691, 0.814373, 0.890479,
    0.422295, 0.821101, 0.894376,
    0.439899, 0.827829, 0.898274,
    0.457504, 0.834557, 0.902172,
    0.475108, 0.841285, 0.90607,
    0.492712, 0.848013, 0.909967,
    0.510316, 0.85474, 0.913865,
    0.527921, 0.861468, 0.917763,
    0.545525, 0.868196, 0.92166,
    0.563129, 0.874924, 0.925558,
    0.580733, 0.881652, 0.929456,
    0.598338, 0.88838, 0.933354,
    0.615738, 0.892714, 0.932322,
    0.632876, 0.89397, 0.924953,
    0.650014, 0.895226, 0.917583,
    0.667152, 0.896482, 0.910214,
    0.68429, 0.897738, 0.902844,
    0.701428, 0.898994, 0.895475,
    0.718566, 0.900251, 0.888106,
    0.735704, 0.901507, 0.880736,
    0.752842, 0.902763, 0.873367,
    0.769979, 0.904019, 0.865998,
    0.787117, 0.905275, 0.858628,
    0.804255, 0.906531, 0.851259,
    0.821393, 0.907787, 0.84389,
    0.838531, 0.909044, 0.83652,
    0.855669, 0.9103, 0.829151,
    0.872807, 0.911556, 0.821782,
    0.883573, 0.909782, 0.806526,
    0.887967, 0.904977, 0.783384,
    0.89236, 0.900173, 0.760242,
    0.896754, 0.895369, 0.7371,
    0.901147, 0.890565, 0.713958,
    0.905541, 0.88576, 0.690816,
    0.909935, 0.880956, 0.667674,
    0.914328, 0.876151, 0.644532,
    0.918722, 0.871347, 0.62139,
    0.923116, 0.866543, 0.598248,
    0.92751, 0.861738, 0.575106,
    0.931903, 0.856934, 0.551964,
    0.936297, 0.85213, 0.528822,
    0.940691, 0.847326, 0.505679,
    0.945084, 0.842521, 0.482538,
    0.949478, 0.837717, 0.459395,
    0.949744, 0.830493, 0.433414,
    0.9468, 0.821387, 0.405225,
    0.943856, 0.81228, 0.377035,
    0.940912, 0.803174, 0.348846,
    0.937968, 0.794068, 0.320656,
    0.935024, 0.784962, 0.292466,
    0.93208, 0.775856, 0.264277,
    0.929136, 0.76675, 0.236087,
    0.926191, 0.757644, 0.207898,
    0.923247, 0.748538, 0.179708,
    0.920303, 0.739431, 0.151519,
    0.917359, 0.730325, 0.123329,
    0.914415, 0.721219, 0.0951396,
    0.911471, 0.712113, 0.06695,
    0.908527, 0.703007, 0.0387605,
    0.905583, 0.693901, 0.0105709,
    0.902502, 0.683442, 0,
    0.899339, 0.672171, 0,
    0.896175, 0.6609, 0,
    0.893012, 0.649629, 0,
    0.889848, 0.638358, 0,
    0.886685, 0.627087, 0,
    0.883522, 0.615817, 0,
    0.880358, 0.604546, 0,
    0.877195, 0.593275, 0,
    0.874032, 0.582004, 0,
    0.870868, 0.570733, 0,
    0.867705, 0.559462, 0,
    0.864542, 0.548192, 0,
    0.861378, 0.536921, 0,
    0.858215, 0.52565, 0,
    0.855051, 0.514379, 0,
    0.850743, 0.503063, 3.81777e-05,
    0.845914, 0.491726, 9.37086e-05,
    0.841084, 0.480389, 0.000149239,
    0.836255, 0.469052, 0.00020477,
    0.831426, 0.457715, 0.000260301,
    0.826596, 0.446378, 0.000315832,
    0.821767, 0.435041, 0.000371363,
    0.816938, 0.423704, 0.000426894,
    0.812108, 0.412368, 0.000482425,
    0.807279, 0.401031, 0.000537955,
    0.80245, 0.389694, 0.000593486,
    0.797621, 0.378357, 0.000649017,
    0.792791, 0.36702, 0.000704548,
    0.787962, 0.355683, 0.000760079,
    0.783133, 0.344346, 0.00081561,
    0.778303, 0.333009, 0.00087114,
    0.772191, 0.321182, 0.000970476,
    0.765651, 0.309192, 0.00108441,
    0.75911, 0.297203, 0.00119834,
    0.75257, 0.285213, 0.00131228,
    0.74603, 0.273223, 0.00142621,
    0.73949, 0.261233, 0.00154015,
    0.73295, 0.249243, 0.00165408,
    0.72641, 0.237253, 0.00176802,
    0.719869, 0.225263, 0.00188195,
    0.713329, 0.213273, 0.00199588,
    0.706789, 0.201283, 0.00210982,
    0.700249, 0.189293, 0.00222375,
    0.693709, 0.177304, 0.00233769,
    0.687168, 0.165314, 0.00245162,
    0.680628, 0.153324, 0.00256556,
    0.674088, 0.141334, 0.00267949,
    0.664498, 0.131995, 0.00256316,
    0.654205, 0.123268, 0.00239369,
    0.643912, 0.114541, 0.00222423,
    0.633618, 0.105814, 0.00205476,
    0.623325, 0.0970873, 0.0018853,
    0.613032, 0.0883604, 0.00171583,
    0.602738, 0.0796334, 0.00154637,
    0.592445, 0.0709064, 0.0013769,
    0.582152, 0.0621795, 0.00120744,
    0.571859, 0.0534525, 0.00103797,
    0.561565, 0.0447255, 0.000868506,
    0.551272, 0.0359986, 0.000699041,
    0.540978, 0.0272716, 0.000529576,
    0.530685, 0.0185446, 0.00036011,
    0.520392, 0.00981769, 0.000190645,
    0.510099, 0.00109072, 2.11802e-05,
    0.497315, 2.01064e-05, 3.01596e-05,
    0.484177, 4.30847e-05, 6.46271e-05,
    0.471038, 6.60631e-05, 9.90945e-05,
    0.457899, 8.90414e-05, 0.000133562,
    0.44476, 0.00011202, 0.000168029,
    0.431622, 0.000134998, 0.000202497,
    0.418483, 0.000157976, 0.000236964,
    0.405344, 0.000180955, 0.000271432,
    0.392205, 0.000203933, 0.000305899,
    0.379067, 0.000226911, 0.000340367,
    0.365928, 0.00024989, 0.000374834,
    0.352789, 0.000272868, 0.000409302,
    0.33965, 0.000295846, 0.000443769,
    0.326512, 0.000318825, 0.000478236,
    0.313373, 0.000341803, 0.000512704,
    0.300234, 0.000364781, 0.000547171,
    0.282726, 0.000540353, 0.000517011,
    0.264928, 0.000726094, 0.000482544,
    0.247129, 0.000911835, 0.000448076,
    0.229331, 0.00109758, 0.000413609,
    0.211532, 0.00128332, 0.000379141,
    0.193733, 0.00146906, 0.000344674,
    0.175935, 0.0016548, 0.000310207,
    0.158136, 0.00184054, 0.000275739,
    0.140337, 0.00202628, 0.000241272,
    0.122539, 0.00221202, 0.000206804,
    0.10474, 0.00239777, 0.000172337,
    0.0869416, 0.00258351, 0.000137869,
    0.0691429, 0.00276925, 0.000103402,
    0.0513443, 0.00295499, 6.89344e-05,
    0.0335457, 0.00314073, 3.44669e-05,
    0.0157473, 0.00332647, 0,
};

static std::vector<float> opacities = {
    0.108614, 
    0.152417, 
    0.19622, 
    0.240023, 
    0.283826, 
    0.327629, 
    0.371432, 
    0.415235, 
    0.459038, 
    0.502841, 
    0.546644, 
    0.590447, 
    0.63425, 
    0.678053, 
    0.721856, 
    0.765659, 
    0.809462, 
    0.853265, 
    0.897068, 
    0.940871, 
    0.957192, 
    0.935508, 
    0.913824, 
    0.89214, 
    0.870456, 
    0.848772, 
    0.827088, 
    0.805404, 
    0.78372, 
    0.762035, 
    0.740351, 
    0.718667, 
    0.696983, 
    0.675299, 
    0.653615, 
    0.631931, 
    0.610247, 
    0.588562, 
    0.566878, 
    0.545194, 
    0.52351, 
    0.501826, 
    0.480142, 
    0.458458, 
    0.436774, 
    0.415089, 
    0.393405, 
    0.371721, 
    0.350037, 
    0.328353, 
    0.306669, 
    0.284985, 
    0.2633, 
    0.241616, 
    0.219932, 
    0.198248, 
    0.176564, 
    0.15488, 
    0.133196, 
    0.111512, 
    0.0725561, 
    0.00127251, 
    0.00807819, 
    0.0163032, 
    0.0245282, 
    0.0327532, 
    0.0409782, 
    0.0492033, 
    0.0574283, 
    0.0656533, 
    0.0738783, 
    0.0821033, 
    0.0903283, 
    0.0985534, 
    0.106778, 
    0.115003, 
    0.123228, 
    0.131453, 
    0.146022, 
    0.201612, 
    0.257202, 
    0.312791, 
    0.368381, 
    0.423971, 
    0.479561, 
    0.53515, 
    0.59074, 
    0.64633, 
    0.701919, 
    0.757509, 
    0.813099, 
    0.868689, 
    0.87781, 
    0.883357, 
    0.888904, 
    0.894451, 
    0.899999, 
    0.905546, 
    0.911093, 
    0.91664, 
    0.922187, 
    0.927734, 
    0.933281, 
    0.938828, 
    0.944375, 
    0.949922, 
    0.95547, 
    0.961017, 
    0.966564, 
    0.972111, 
    0.977658, 
    0.98156, 
    0.982382, 
    0.983205, 
    0.984027, 
    0.98485, 
    0.985672, 
    0.986495, 
    0.987317, 
    0.98814, 
    0.988962, 
    0.989785, 
    0.990607, 
    0.99143, 
    0.992252, 
    0.993075, 
    0.993897, 
    0.99472, 
    0.995542, 
    0.996365, 
    0.997187, 
    0.99801, 
    0.998832, 
    0.999655, 
    0.99947, 
    0.998556, 
    0.997642, 
    0.996728, 
    0.99375, 
    0.988555, 
    0.983361, 
    0.978166, 
    0.972971, 
    0.967776, 
    0.962582, 
    0.957387, 
    0.952192, 
    0.980906, 
    0.999888, 
    0.999646, 
    0.999404, 
    0.999162, 
    0.99892, 
    0.998678, 
    0.998436, 
    0.998194, 
    0.997952, 
    0.99771, 
    0.997469, 
    0.997227, 
    0.996985, 
    0.996743, 
    0.996501, 
    0.996259, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0.00290816, 
    0.0149715, 
    0.0270348, 
    0.0390982, 
    0.0511615, 
    0.0632249, 
    0.0752882, 
    0.187601, 
    0.443752, 
    0.699902, 
    0.67856, 
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
    0.0176842, 
    0.0412626, 
    0.0648409, 
    0.0884193, 
    0.111998, 
    0.135576, 
    0.159154, 
    0.139365, 
    0.115786, 
    0.0922078, 
    0.0686294, 
    0.045051, 
    0.0214726, 
    7.34577e-05, 
    0.000895956, 
    0.00171846, 
    0.00254095, 
    0.00336345, 
    0.127564, 
    0.358687, 
    0, 
    0, 
    0, 
    0, 
    0, 
    0.00645445, 
    0.0156048, 
    0.0247551, 
    0.0339054, 
    0.0430558, 
    0.0522061, 
    0.0613564, 
    0.0610811, 
    0.0477716, 
    0.034462, 
    0.0211525, 
    0.00784295, 
    0.00144782, 
    0.00497283, 
    0.00849783, 
    0.00572768, 
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


OSPGeometry make_triangle() {
  // Data for a triangle
  // const float tri_verts[] = {-4, -4, 0, 0, -4, 4, 0, 0, 4, -4, 0, 0};
  const float tri_verts[] = {11.995, 12, 0.1, 0, 11.995, 20, 0.1, 0, 19.995, 12, 0.1, 0};

	const float tri_color[] =  {
		1, 0, 0, 1,
		0, 1, 0, 1,
		0, 0, 1, 1
	};
	const int32_t tri_index[] = { 0, 1, 2 };

	// create and setup model and mesh
	OSPGeometry mesh = ospNewGeometry("triangles");
	OSPData tri_data = ospNewData(3, OSP_FLOAT3A, tri_verts, 0);
	ospCommit(tri_data);
	ospSetData(mesh, "vertex", tri_data);
	ospRelease(tri_data);

	tri_data = ospNewData(3, OSP_FLOAT4, tri_color, 0);
	ospCommit(tri_data);
	ospSetData(mesh, "vertex.color", tri_data);
	ospRelease(tri_data);

	tri_data = ospNewData(1, OSP_INT3, tri_index, 0);
	ospCommit(tri_data);
	ospSetData(mesh, "index", tri_data);
	ospRelease(tri_data);
	//ospSetMaterial(mesh, mat);
	ospCommit(mesh);
	return mesh;
}


std::array<osp::vec4f, 3>  isoColors = {osp::vec4f{254.0f/255.f,129.0/255.f,0.0/4.f,1.0f},
osp::vec4f{0.0f/255.f,98.0/255.f,254.0/255.f,0.6f},osp::vec4f{0.0f,0.0f,1.0f,0.6f}};

static ospcommon::vec2f valueRange{0.f, -1.f};

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
  if (ospInit(&ac, av) != OSP_NO_ERROR) {
    throw std::runtime_error("FATAL ERROR DURING INITIALIZATION!");
    return 1;
  }

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
    else if (str == "-valueRange") {
      ospray::impi::Parse<2>(ac, av, i, valueRange);
    }
    else if (str == "-iso" || str == "-isoValue") {
      ospray::impi::Parse<1>(ac, av, i, isoValues[0].v);
    }
    else if (str == "-isos" || str == "-isoValues") {
      try {
	int n = 0;
	ospray::impi::Parse<1>(ac, av, i, n);
	isoValues.resize(n);
	for (int j = 0; j < n; ++j) {
	  ospray::impi::Parse<1>(ac, av, i, isoValues[j].v);
	}
      } catch (const std::runtime_error& e) {
	throw std::runtime_error(std::string(e.what())+
				 " usage: -isos "
				 "<# of iso-values> "
				 "<iso-value list>");
      }
    }
    else if (str == "-isoColor") {
      ospray::impi::Parse<3>(ac, av, i, isoValues[0].c);
    }
    else if (str == "-isoColors") {
      try {
	for (int j = 0; j < isoValues.size(); ++j) {
	  ospray::impi::Parse<3>(ac, av, i, isoValues[j].c);
	}
      } catch (const std::runtime_error& e) {
	throw std::runtime_error(std::string(e.what())+
				 " usage: -isoColors "
				 "<R, G, B list>");
      }
    }
    else if (str == "-translate") {
      ospray::impi::Parse<3>(ac, av, i, objTranslate);
    }
    else if (str == "-scale") {
      ospray::impi::Parse<3>(ac, av, i, objScale);
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
  auto device = ospGetCurrentDevice();
  if (device == nullptr) {
    throw std::runtime_error("FATAL ERROR DURING GETTING CURRENT DEVICE!");
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
  if (ospLoadModule("impi") != OSP_NO_ERROR) {
    throw std::runtime_error("failed to initialize IMPI module");
  }
  
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
  if (valueRange.x > valueRange.y) {
    ospSetVec2f(transferFcn, "valueRange", 
		(const osp::vec2f&)amrVolume->Range());
  } else {
    ospSetVec2f(transferFcn, "valueRange", 
		(const osp::vec2f&)valueRange);
  }
  ospCommit(transferFcn);
  ospRelease(colorsData);
  ospRelease(opacitiesData);

  // setup volume
  OSPVolume volume = amrVolume->Create(transferFcn);
  if (showVolume) {
    ospAddVolume(world, volume);
  }

  // setup isosurfaces
  switch (isoMode) {
  case NORMAL:
    // --> normal isosurface
    {
      // Node: all isosurfaces will share one material for now
      OSPMaterial nmtl = ospNewMaterial(renderer, "OBJMaterial");
      ospSetVec3f(nmtl, "Kd", osp::vec3f{1.0f, 1.0f, 1.0f});
      ospSetVec3f(nmtl, "Ks", osp::vec3f{0.1f, 0.1f, 0.1f});
      ospSet1f(nmtl, "Ns", 10.f);
      ospCommit(nmtl);
      std::vector<float> vlist(isoValues.size());
      for (auto i = 0; i < isoValues.size(); ++i) {
	vlist[i] = isoValues[i].v;
      }
      OSPGeometry niso = ospNewGeometry("isosurfaces");
      OSPData niso_values = ospNewData(vlist.size(), OSP_FLOAT, 
				       vlist.data());
      ospSetData(niso, "isovalues", niso_values);
      ospSetObject(niso, "volume", volume);
      ospSetMaterial(niso, nmtl); // see performance impact (x7 slower for cosmos)
      ospCommit(niso);
      ospAddGeometry(world, niso);
    }
    break;
  case IMPI:
    // --> implicit isosurface
    {
      // Note: because there is no naive multi-iso surface support,
      //       we build multiple iso-geometries here      
      for (auto& v : isoValues) {
	std::cout << "v = " << v.v << " "
		  << "c = " << v.c.x << " " << v.c.y << " " << v.c.z
		  << std::endl;
	if (rendererName == "scivis") {
	  v.mtl = ospNewMaterial(renderer, "OBJMaterial");
	  ospSetVec3f(v.mtl, "Kd", (const osp::vec3f&)v.c);
	  ospSetVec3f(v.mtl, "Ks", osp::vec3f{0.1f, 0.1f, 0.1f});
	  ospSet1f(v.mtl, "Ns", 10.f);
	  ospCommit(v.mtl);
	} else {
	  //----------------
	  //v.mtl = ospNewMaterial(renderer, "ThinGlass");
	  //ospSetVec3f(v.mtl, "attenuationColor", (const osp::vec3f&)v.c);
	  //ospSet1f(v.mtl, "thickness", 0.1f);
	  //----------------
	  //v.mtl = ospNewMaterial(renderer, "Alloy");
	  //ospSetVec3f(v.mtl, "color", (const osp::vec3f&)v.c);
	  //----------------
	  v.mtl = ospNewMaterial(renderer, "MetallicPaint");
	  ospSetVec3f(v.mtl, "baseColor", (const osp::vec3f&)v.c);
	  //----------------
	  ospCommit(v.mtl);
	}
	v.geo = ospNewGeometry("impi"); 
	ospSet1f(v.geo, "isoValue", v.v);
	ospSetObject(v.geo, "amrDataPtr", volume);
	ospSetMaterial(v.geo, v.mtl); // see performance impact (x7 slower for cosmos)
	ospCommit(v.geo);
	ospAddGeometry(world, v.geo);
      }
    }
    break;
  default:
    throw std::runtime_error("wrong ISO-Mode, this shouldn't happen");
  }
  
  // debugging
  //OSPGeometry triangle = make_triangle();
  //ospAddGeometry(world, triangle);

  // setup object
  Mesh mesh;
  affine3f transform = 
    affine3f::translate(objTranslate) * 
    affine3f::scale(objScale);
  if (showObject) {
    OSPMaterial mtlobj = ospNewMaterial(renderer, "OBJMaterial");
    ospSetVec3f(mtlobj, "Kd", osp::vec3f{3/255.f, 10/255.f, 25/255.f});
    ospSetVec3f(mtlobj, "Ks", osp::vec3f{77/255.f, 77/255.f, 77/255.f});
    ospSet1f(mtlobj, "Ns", 10.f);
    ospCommit(mtlobj);
    mesh.LoadFromFileObj(inputMesh.c_str(), false);
    mesh.SetTransform(transform);
    mesh.AddToModel(world, renderer, mtlobj);
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
	      osp::vec3f{0.f, 0.f, 0.f});
  ospSetData(renderer, "lights", lights);
  ospSetObject(renderer, "model", world);
  ospSetObject(renderer, "camera", camera);
  ospSet1i(renderer, "shadowEnabled", 1);
  ospSet1i(renderer, "oneSidedLighting", 1);
  ospSet1i(renderer, "maxDepth", 100);
  ospSet1i(renderer, "spp", 1);
  ospSet1i(renderer, "autoEpsilon", 1);
  ospSet1i(renderer, "aoSamples", 1);
  ospSet1i(renderer, "aoTransparencyEnabled", 1);
  ospSet1f(renderer, "aoDistance", 10000.0f);
  ospSet1f(renderer, "epsilon", 0.001f);
  ospSet1f(renderer, "minContribution", 0.001f);
  ospCommit(renderer);

#if USE_VIEWER

  ospray::viewer::Handler(camera, 
			  (const osp::vec3f&)vp, 
			  (const osp::vec3f&)vu, 
			  (const osp::vec3f&)vi);
  ospray::viewer::Handler(transferFcn, amrVolume->Range().x, amrVolume->Range().y);
  ospray::viewer::Handler(world, renderer);
  ospray::viewer::Render(window);

#else

  // setup framebuffer
  OSPFrameBuffer fb = ospNewFrameBuffer((const osp::vec2i&)imgSize, 
					OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
  ospFrameBufferClear(fb, OSP_FB_COLOR | OSP_FB_ACCUM);

  // render 10 more frames, which are accumulated to result in a better converged image
  std::cout << "#osp:bench: start warmups for " 
	    << numFrames.x << " frames" << std::endl;
  for (int frames = 0; frames < numFrames.x; frames++) { // skip some frames to warmup
    ospRenderFrame(fb, renderer, OSP_FB_COLOR | OSP_FB_ACCUM);
  }
  std::cout << "#osp:bench: done warmups" << std::endl;
  std::cout << "#osp:bench: start benchmarking for "
	    << numFrames.y << " frames" << std::endl;
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
