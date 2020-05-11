#include "gui_windows.h"


bool ARTSGUI::Windows::full(GLFWwindow* window, const ImVec2 origpos, const char* name)
{
  return sub<1, 1, 0, 0, 1, 1>(window, origpos, name);
}
