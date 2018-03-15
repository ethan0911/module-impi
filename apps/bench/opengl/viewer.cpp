// ======================================================================== //
// Copyright SCI Institute, University of Utah, 2018
// ======================================================================== //
#include "viewer.h"
#include "common.h"
#include "camera.h"
#include "framebuffer.h"

#include <thread>
#include <atomic>

#include <imgui.h>
#include <imgui_glfw_impi.h>

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
static std::vector<GLFWwindow*> windowmap;
static Camera      camera;
static Framebuffer framebuffer;
static OSPCamera           ospCam;
static OSPTransferFunction ospTfn;
static OSPModel            ospMod;
static OSPRenderer         ospRen;
static bool printTFN = false;

static std::vector<float> cptr;
static std::vector<float> aptr;
static bool tfnptr = false;

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
static affine3f Identity(vec3f(1,0,0), vec3f(0,1,0), vec3f(0,0,1), vec3f(0,0,0));
struct Sphere {
  struct SphereInfo
  {
    vec3f org{0.f,0.f,0.f};
    int   colorID{0};
  };
  OSPGeometry  sphere;
  OSPGeometry  instance;
  OSPData      data;
  OSPModel     local;
  SphereInfo   info;
  bool added = false;
  void Init()
  {
    data = ospNewData(sizeof(SphereInfo), OSP_UCHAR, &info, OSP_DATA_SHARED_BUFFER);
    ospCommit(data);
    vec4f color(1.f, 0.f, 0.f, 1.f);
    OSPData cdata = ospNewData(1, OSP_FLOAT4, &color);
    ospCommit(cdata);
    sphere = ospNewGeometry("spheres");
    ospSetData(sphere, "spheres", data);
    ospSetData(sphere, "color", cdata);
    ospSet1i(sphere, "offset_center", 0);
    ospSet1i(sphere, "offset_colorID", int(sizeof(vec3f)));
    ospSet1f(sphere, "radius", 0.01f * camera.CameraFocalLength());
    ospCommit(sphere);
    ospRelease(cdata);
    local = ospNewModel();
    ospAddGeometry(local, sphere);
    ospCommit(local);
    instance = ospNewInstance(local, (osp::affine3f &)Identity);
    ospCommit(instance);
  }
  void Update(const vec3f& center, OSPModel mod)
  {
    if (added) {
      if (info.org != center) {
	info.org = center;
	ospCommit(sphere);
	ospCommit(local);
	ospCommit(mod);
      }
    }
  }
  void Add(OSPModel mod)
  {
    added = true;
    ospAddGeometry(mod, instance);
    ospCommit(mod);
  }
  void Remove(OSPModel mod)
  {
    ospRemoveGeometry(mod, instance);
    ospCommit(mod);
    added = false;
  }
} sphere;

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
#include "widgets/TransferFunctionWidget.h"
static std::shared_ptr<tfn::tfn_widget::TransferFunctionWidget> tfnWidget;
static float tfnValueRange[2] = {0.f, 1.f};

namespace ospray {

  void RenderWindow(GLFWwindow *window);
  GLFWwindow *CreateWindow();

  namespace viewer {
    int Init(const int ac, const char** av, const size_t& w, const size_t& h) 
    {
      camera.SetSize(w, h);
      windowmap.push_back(CreateWindow());
      return windowmap.size() - 1;
    }
    void Render(int id) { 
      sphere.Init();
      framebuffer.Init(camera.CameraWidth(), camera.CameraHeight(), ospRen);
      RenderWindow(windowmap[id]); 
    };
    void Handler(OSPCamera c, const osp::vec3f& vp, const osp::vec3f& vu, const osp::vec3f& vi)
    { 
      ospCam = c; 
      camera.SetViewPort(vec3f(vp.x, vp.y, vp.z), 
			 vec3f(vu.x, vu.y, vu.z), 
			 vec3f(vi.x, vi.y, vi.z), 60.f);
      camera.Init(ospCam);
    };
    void Handler(OSPTransferFunction t, const float& a, const float& b) 
    { 
      ospTfn = t; 
      tfnValueRange[0] = a;
      tfnValueRange[1] = b;
    };
    void Handler(OSPModel m, OSPRenderer r) { ospMod = m; ospRen = r; };
  };

  //-----------------------------------------------------------------------------------------//
  static std::thread *osprayThread = nullptr;
  static std::atomic<bool> osprayStop(false);
  static std::atomic<bool> osprayClear(false);
  void StartOSPRay() {
    osprayStop = false;
    osprayThread = new std::thread([=] {
	while (!osprayStop) {
	  if (osprayClear) { framebuffer.Clear(); osprayClear = false; }
	  sphere.Update(camera.CameraFocus(), ospMod);
	  framebuffer.Render();
	}
      });
  }
  void StopOSPRay() { osprayStop = true; osprayThread->join(); delete osprayThread; }
  void ClearOSPRay() { osprayClear = true; }
  void ResizeOSPRay(int width, int height)
  {
    StopOSPRay();
    framebuffer.Resize((size_t)width, (size_t)height);
    StartOSPRay();
  }
  void UploadOSPRay()
  {
    framebuffer.Display();
  }

  //-----------------------------------------------------------------------------------------//

  void error_callback(int error, const char *description) {
    fprintf(stderr, "Error: %s\n", description);
  }

  void char_callback(GLFWwindow *window, unsigned int c) {
    ImGuiIO& io = ImGui::GetIO();
    if (c > 0 && c < 0x10000) { io.AddInputCharacter((unsigned short)c); }
  }

  void key_onhold_callback(GLFWwindow *window) {

    if (ImGui::GetIO().WantCaptureKeyboard) {
      return;
    }

    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {  
      /* UP: forward */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMoveNZ(1.f * camera.CameraFocalLength());
      } else {
	camera.CameraMoveNZ(0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();
    }
    else if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
      /* DOWN: backward */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMoveNZ(-0.5f * camera.CameraFocalLength());
      } else {
	camera.CameraMoveNZ(-0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();      
    }
    else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
      /* A: left */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMovePX(1.f * camera.CameraFocalLength());
      } else {
	camera.CameraMovePX(0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();
    }
    else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
      /* D: right */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMovePX(-0.5f * camera.CameraFocalLength());
      } else {
	camera.CameraMovePX(-0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();      
    }
    else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
      /* S: down */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMovePY(1.f * camera.CameraFocalLength());
      } else {
	camera.CameraMovePY(0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();
    }
    else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
      /* W: up */
      if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
	camera.CameraMovePY(-0.5f * camera.CameraFocalLength());
      } else {
	camera.CameraMovePY(-0.01f * camera.CameraFocalLength());
      }
      ClearOSPRay();      
    }

  }

  void key_onpress_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, GL_TRUE);
    }     
    if (!ImGui::GetIO().WantCaptureKeyboard) {
      if (key == GLFW_KEY_LEFT_ALT) {
	StopOSPRay();
	if (action == GLFW_PRESS) {
	  sphere.Add(ospMod);
	} else if (action == GLFW_RELEASE) {
	  sphere.Remove(ospMod);
	}
	ClearOSPRay();
	StartOSPRay();
      }
      if (key == GLFW_KEY_P && action == GLFW_PRESS) {
	if (tfnptr) {
	  const std::vector<float> & c = cptr;
	  const std::vector<float> & a = aptr;
	  std::cout << std::endl << "static std::vector<float> colors = {" << std::endl;
	  for (int i = 0; i < c.size()/3; ++i) {
	    std::cout << "    " << c[3 * i] << ", " << c[3 * i + 1] << ", " << c[3 * i + 2] << "," << std::endl;
	  }
	  std::cout << "};" << std::endl;
	  std::cout << "static std::vector<float> opacities = {" << std::endl;
	  for (int i = 0; i < a.size()/2; ++i) {
	    std::cout << "    " << a[2 * i + 1] << ", " << std::endl;
	  }
	  std::cout << "};" << std::endl << std::endl;
	}
      }
      else if (key == GLFW_KEY_V && action == GLFW_PRESS) {
	const auto vi = camera.CameraFocus();
	const auto vp = camera.CameraPos();
	const auto vu = camera.CameraUp();
	std::cout << std::endl 
		  << "-vp " << vp.x << " " << vp.y << " " << vp.z << " "
		  << "-vi " << vi.x << " " << vi.y << " " << vi.z << " "
		  << "-vu " << vu.x << " " << vu.y << " " << vu.z << " "
		  << std::endl;
      }
    } else {      
      ImGui_Impi_KeyCallback(window, key, scancode, action, mods); 
    }
  }

  void cursor_position_callback(GLFWwindow *window, double xpos, double ypos) {
    if (!ImGui::GetIO().WantCaptureMouse) {
      int left_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
      int right_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
      if (left_state == GLFW_PRESS) {
	camera.CameraDrag((float)xpos, (float)ypos);
	ClearOSPRay();
      } else {
	camera.CameraBeginDrag((float)xpos, (float)ypos);
      }
      if (right_state == GLFW_PRESS) {
	camera.CameraZoom((float)xpos, (float)ypos);
	ClearOSPRay();
      } else {
	camera.CameraBeginZoom((float)xpos, (float)ypos);
      }
    }
  }

  void window_size_callback(GLFWwindow *window, int width, int height) {
    glViewport(0, 0, width, height);
    camera.CameraUpdateProj((size_t)width, (size_t)height);
    ResizeOSPRay(width, height);
  }

  //-----------------------------------------------------------------------------------------//

  void RenderWindow(GLFWwindow *window)
  {
    // Init   
    ImGui_Impi_Init(window, false);
    if (ospTfn != nullptr) {
      tfnWidget = std::make_shared<tfn::tfn_widget::TransferFunctionWidget>
	([&](const std::vector<float> &c, 
	     const std::vector<float> &a,
	     const std::array<float, 2>& r) 
	 {
	   tfnptr = true;	
	   cptr = std::vector<float>(c);
	   aptr = std::vector<float>(a);	 
	   OSPData colorsData = ospNewData(c.size() / 3, OSP_FLOAT3, c.data());
	   ospCommit(colorsData);
	   std::vector<float>o(a.size()/2);
	   for (int i = 0; i < a.size()/2; ++i) { o[i] = a[2*i+1]; }
	   OSPData opacitiesData = ospNewData(o.size(), OSP_FLOAT, o.data());
	   ospCommit(opacitiesData);
	   ospSetData(ospTfn, "colors",    colorsData);
	   ospSetData(ospTfn, "opacities", opacitiesData);
	   ospSetVec2f(ospTfn, "valueRange", osp::vec2f{r[0], r[1]});
	   ospCommit(ospTfn);
	   ospRelease(colorsData);
	   ospRelease(opacitiesData);
	   ClearOSPRay();
	 });
      tfnWidget->setDefaultRange(tfnValueRange[0], tfnValueRange[1]);
    }
    // Start
    StartOSPRay();
    while (!glfwWindowShouldClose(window)) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      {
	key_onhold_callback(window);
	UploadOSPRay();
	ImGui_Impi_NewFrame();
	if (ospTfn != nullptr) {
	  if (tfnWidget->drawUI()) { tfnWidget->render(); };
	}
	ImGui::Render();
      }
      glfwSwapBuffers(window);
      glfwPollEvents();
    }
    // ShutDown
    StopOSPRay();
    {
      ImGui_Impi_Shutdown();
    }
    glfwDestroyWindow(window);
    glfwTerminate();
  }

  //-----------------------------------------------------------------------------------------//

  GLFWwindow *CreateWindow() {
    // Initialize GLFW
    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) { exit(EXIT_FAILURE); }
    // Provide Window Hints
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Create Window
    GLFWwindow *window = glfwCreateWindow(camera.CameraWidth(),
					  camera.CameraHeight(),
					  "OSPRay Test Viewer",
					  nullptr, nullptr);
    if (!window) {
      glfwTerminate();
      exit(EXIT_FAILURE);
    }
    // Callback
    glfwSetKeyCallback(window, key_onpress_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetCharCallback(window, char_callback);
    // Ready
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
    glfwSwapInterval(1);
    check_error_gl("Ready");
    // Setup OpenGL
    glEnable(GL_DEPTH_TEST);
    check_error_gl("Setup OpenGL Options");
    return window;
  }

};
