#ifndef GUI_WINDOWS_H
#define GUI_WINDOWS_H

#include "gui_inc.h"

namespace ARTSGUI {
namespace Windows {

template <unsigned WIDTH, unsigned HEIGHT, unsigned WIDTH_POS, unsigned HEIGHT_POS>
bool sub(GLFWwindow* window, const ImVec2 origpos, const char * name)
{
  static_assert(WIDTH, "No size window not allowed");
  static_assert(HEIGHT, "No size window not allowed");
  static_assert(WIDTH > WIDTH_POS, "More width than possible");
  static_assert(HEIGHT > HEIGHT_POS, "More height than possible");
  constexpr float wscale = 1.0f / float(WIDTH);
  constexpr float hscale = 1.0f / float(HEIGHT);
  
  //Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 size = {float(width)*wscale-2*origpos.x, float(height)*hscale-origpos.y};
  ImVec2 pos = {origpos.x + size.x*WIDTH_POS, origpos.y + size.y*HEIGHT_POS};
  
  // Show a simple window
  ImGui::SetNextWindowPos(pos);
  ImGui::SetNextWindowSize(size);
  
  return ImGui::Begin(name, NULL, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar);
}

bool full(GLFWwindow* window, const ImVec2 origpos, const char * name);

};  // Windows
};  // ARTSGUI

#endif  // GUI_WINDOWS_H
