#ifndef GUI_MENUBAR_H
#define GUI_MENUBAR_H

#include "gui_inc.h"

namespace ARTSGUI {
namespace MainMenu {
struct Config {
  ImGuiIO& io;
  bool show_about_help, show_metrics_help, show_style_help;
  bool fullscreen;
  int width, height, xpos, ypos;
  Config(bool fullscreen_on=false) : io(ImGui::GetIO()),
  show_about_help(false), show_metrics_help(false), show_style_help(false), 
  fullscreen(fullscreen_on), width(1280), height(720), xpos(50), ypos(50)
  {
    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
  }
};  // Config

void fullscreen(Config& cfg, GLFWwindow* window);
void quitscreen(Config& cfg, GLFWwindow* window);
void imgui_help(Config& cfg);
void arts_help();

};  // MainMenu
};  // ARTSGUI

#endif  // GUI_MENUBAR_H
