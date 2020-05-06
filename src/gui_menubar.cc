#include "gui_menubar.h"
#include "gui_help.h"
#include "auto_version.h"

void ARTSGUI::MainMenu::fullscreen(Config& cfg, GLFWwindow* window)
{
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("Fullscreen", "F11")) {
        if (not cfg.fullscreen) {
          glfwGetWindowSize(window, &cfg.width, &cfg.height);
          glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
          const GLFWvidmode * mode = glfwGetVideoMode(get_current_monitor(window));
          glfwSetWindowMonitor(window, get_current_monitor(window), 0, 0, mode->width, mode->height, 0);
        }
        else
          glfwSetWindowMonitor(window, NULL, cfg.xpos, cfg.ypos, cfg.width, cfg.height, 0);
        
        cfg.fullscreen = not cfg.fullscreen;
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (ImGui::IsKeyPressed(GLFW_KEY_F11) or (cfg.fullscreen and ImGui::IsKeyPressed(GLFW_KEY_ESCAPE))) {
    if (not cfg.fullscreen) {
      glfwGetWindowSize(window, &cfg.width, &cfg.height);
      glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
      const GLFWvidmode * mode = glfwGetVideoMode(get_current_monitor(window));
      glfwSetWindowMonitor(window, get_current_monitor(window), 0, 0, mode->width, mode->height, 0);
    }
    else
      glfwSetWindowMonitor(window, NULL, cfg.xpos, cfg.ypos, cfg.width, cfg.height, 0);
    cfg.fullscreen = not cfg.fullscreen;
  }
}

void ARTSGUI::MainMenu::quitscreen(Config& cfg, GLFWwindow* window)
{
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("Quit", "Ctrl+X"))
        glfwSetWindowShouldClose(window, 1);
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (cfg.io.KeyCtrl and ImGui::IsKeyPressed(GLFW_KEY_X)) {
    glfwSetWindowShouldClose(window, 1);
  }
}

void ARTSGUI::MainMenu::imgui_help(Config& cfg)
{
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Help")) {
      ImGui::Checkbox("  About ImGui", &cfg.show_about_help);
      ImGui::Checkbox("  Metrics ImGui", &cfg.show_metrics_help);
      ImGui::Checkbox("  Style ImGui", &cfg.show_style_help);
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (cfg.show_about_help) {
    ImGui::ShowAboutWindow();
  }
  if (cfg.show_metrics_help) {
    ImGui::ShowMetricsWindow();
  }
  if (cfg.show_style_help) {
    ImGui::ShowStyleEditor();
  }
}

void ARTSGUI::MainMenu::arts_help()
{
  bool arts_help_popup=false;
  bool arts_license_popup=false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Help")) {
      if (ImGui::MenuItem("About ARTS")) arts_help_popup = true;
      if (ImGui::MenuItem("License ARTS")) arts_license_popup = true;
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (arts_help_popup) {
    ImGui::OpenPopup("ARTS Help");
  }
  
  if (arts_license_popup) {
    ImGui::OpenPopup("ARTS License");
  }
  
  if (ImGui::BeginPopup("ARTS Help")) {
    ImGui::Text("You are using " ARTS_FULL_VERSION);
    ImGui::EndPopup();
  }
  
  if (ImGui::BeginPopup("ARTS License")) {
    ImGui::Text("The ARTS program is free software; you can redistribute it\n"
                "and/or modify it under the terms of the GNU General Public\n"
                "License as published by the Free Software Foundation; either\n"
                "version 2, or (at your option) any later version.\n"
                "\n"
                "This program is distributed in the hope that it will be\n"
                "useful, but WITHOUT ANY WARRANTY; without even the implied\n"
                "warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR\n"
                "PURPOSE. See the GNU General Public License for more\n"
                "details. \n"
                "\n"
                "You should have received a copy of the GNU General Public\n"
                "License along with the program; if not, write to the Free\n"
                "Software Foundation, Inc., 59 Temple Place - Suite 330,\n"
                "Boston, MA 02111-1307, USA.\n");
    ImGui::EndPopup();
  }
}

