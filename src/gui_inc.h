#ifndef GUI_INC_H
#define GUI_INC_H

#include "../3rdparty/gui/imgui/imgui.h"
#include "../3rdparty/gui/imgui/imgui_impl_glfw.h"
#include "../3rdparty/gui/imgui/imgui_impl_opengl3.h"
#include "../3rdparty/gui/implot/implot.h"
#include <stdio.h>

#include <GL/glew.h>     // Initialize with glewInit()
#include <GLFW/glfw3.h>  // Include glfw3.h after our OpenGL definitions

#include "gui_macros.h"

extern "C" {
  inline static void glfw_error_callback(int error, const char* description)
  {
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
  }
}

namespace ARTSGUI {
  void LayoutAndStyleSettings();
};

#endif  // GUI_INC_H
