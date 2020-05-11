#ifndef GUI_WINDOWS_H
#define GUI_WINDOWS_H

#include "gui_inc.h"

namespace ARTSGUI {
namespace Windows {

template <unsigned WIDTH, unsigned HEIGHT, unsigned WIDTH_POS, unsigned HEIGHT_POS, unsigned WIDTH_EXTENT=1, unsigned HEIGHT_EXTENT=1>
bool sub(GLFWwindow* window, const ImVec2 origpos, const char * name)
{
  static_assert(WIDTH, "None size window not allowed");
  static_assert(HEIGHT, "None size window not allowed");
  static_assert(WIDTH_EXTENT, "None extent window not allowed");
  static_assert(HEIGHT_EXTENT, "None extent window not allowed");
  static_assert(WIDTH > WIDTH_POS + WIDTH_EXTENT - 1, "More width than possible");
  static_assert(HEIGHT > HEIGHT_POS + HEIGHT_EXTENT - 1, "More height than possible");
  constexpr float wscale = 1.0f / float(WIDTH);
  constexpr float hscale = 1.0f / float(HEIGHT);
  
  //Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 size = {float(width)*wscale, (float(height)-origpos.y)*hscale};
  ImVec2 pos = {origpos.x + size.x*WIDTH_POS, origpos.y + size.y*HEIGHT_POS};
  size.x *= WIDTH_EXTENT;
  size.y *= HEIGHT_EXTENT;
  
  // Set a simple window frame
  ImGui::SetNextWindowPos(pos);
  ImGui::SetNextWindowSize(size);
  
  return ImGui::Begin(name, NULL, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar);
}

bool full(GLFWwindow* window, const ImVec2 origpos, const char * name);

};  // Windows
};  // ARTSGUI

#endif  // GUI_WINDOWS_H
