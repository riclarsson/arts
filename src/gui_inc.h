#ifndef GUI_INC_H
#define GUI_INC_H

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>
#include <stdio.h>

#include <GL/glew.h>     // Initialize with glewInit()
#include <GLFW/glfw3.h>  // Include glfw3.h after our OpenGL definitions

#include "gui_macros.h"

#include <limits>

extern "C" {
  inline static void glfw_error_callback(int error, const char* description)
  {
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
  }
}

namespace ARTSGUI {
  struct Config {
    ImGuiIO& io;
    bool show_about_help, show_metrics_help, show_style_help;
    bool fullscreen;
    bool autoscale_x;
    int width, height, xpos, ypos;
    float mag_copy;
    Config(bool fullscreen_on=false) : io(ImGui::GetIO()),
    show_about_help(false), show_metrics_help(false), show_style_help(false), 
    fullscreen(fullscreen_on), autoscale_x(false), width(1280), height(720), xpos(50), ypos(50),
    mag_copy(std::numeric_limits<float>::quiet_NaN())
    {
      (void)io;
      io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    }
  };  // Config
  
  void LayoutAndStyleSettings();
};

#endif  // GUI_INC_H
