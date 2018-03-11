#include "viewer.h"
#include "common.h"
#include "camera.h"
#include "framebuffer.h"

#include <thread>
#include <atomic>

#include <imgui.h>
#include <imgui_glfw_impi.h>

#include "widgets/TransferFunctionWidget.h"

static std::shared_ptr<tfn::tfn_widget::TransferFunctionWidget> tfnWidget;
static std::vector<GLFWwindow*> windowmap;

static Camera      camera;
static Framebuffer framebuffer;

static OSPCamera           ospCam;
static OSPTransferFunction ospTfn;
static OSPRenderer         ospRen;

static float tfnValueRange[2] = {0.f, 1.f};

namespace ospray {

  void RenderWindow(GLFWwindow *window);
  GLFWwindow *CreateWindow();

  namespace viewer {
    int Init(const int ac, const char** av) {
      windowmap.push_back(CreateWindow());
      return windowmap.size() - 1;
    }
    void Render(int id) { 
      camera.Init(ospCam, 1.0f);
      framebuffer.Init(camera.CameraWidth(), camera.CameraHeight(), ospRen);
      RenderWindow(windowmap[id]); 
    };
    void Handler(OSPCamera c) { ospCam = c; };
    void Handler(OSPTransferFunction t, const float& a, const float& b) 
    { 
      ospTfn = t; 
      tfnValueRange[0] = a;
      tfnValueRange[1] = b;
    };
    void Handler(OSPRenderer r) { ospRen = r; };
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

  void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, GL_TRUE);
    }
    if (!ImGui::GetIO().WantCaptureKeyboard) {}
    ImGui_Impi_KeyCallback(window, key, scancode, action, mods); 
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
    tfnWidget = std::make_shared<tfn::tfn_widget::TransferFunctionWidget>
      ([&](const std::vector<float> &c, 
	   const std::vector<float> &a,
	   const std::array<float, 2>& r) 
       {
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
    // Start
    StartOSPRay();
    while (!glfwWindowShouldClose(window)) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      {
	UploadOSPRay();
	ImGui_Impi_NewFrame();
	if (tfnWidget->drawUI()) { tfnWidget->render(); };
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
    GLFWwindow *window = glfwCreateWindow((int)camera.CameraWidth(),
					  (int)camera.CameraHeight(),
					  "OSPRay Volume Test Renderer",
					  nullptr, nullptr);
    if (!window) {
      glfwTerminate();
      exit(EXIT_FAILURE);
    }
    // Callback
    glfwSetKeyCallback(window, key_callback);
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
