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


bool ARTSGUI::MainMenu::SelectAtmosphere(VectorView abs_p, VectorView abs_t, MatrixView abs_vmrs, const ArrayOfString& spec_list, const Index level)
{
  bool new_plot=false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Atmospheric properties")) {
      float temperature = float(abs_t[level]);
      if (ImGui::SliderFloat("Temperature [K]", &temperature, 150.0f, 600.0f)) {
        if (temperature < 0.01)
          abs_t[level] = 0.01; 
        else
          abs_t[level] = Numeric(temperature);
        new_plot = true;
      }
      ImGui::Separator();
      
      float pressure = float(abs_p[level]);
      if (ImGui::SliderFloat("Pressure [Pa]", &pressure, 0.01f, 1'000'000.0f, "%g", 10.0f)) {
        if (pressure < 0.01)
          abs_p[level] = 0.01; 
        else
          abs_p[level] = Numeric(pressure);
        new_plot = true;
      }
      ImGui::Separator();
      
      for (Index i=0; i<spec_list.nelem(); i++) {
        float vmr = float(abs_vmrs(i, level) * 1e6);
        if (ImGui::SliderFloat(spec_list[i].c_str(), &vmr, 0.00f, 1'000'000.0f, "%g", 10.0f)) {
          if (vmr < 0)
            abs_vmrs(i, level) = 0;
          else
            abs_vmrs(i, level) = Numeric(vmr) / 1e6;
          new_plot = true;
        }
        ImGui::Separator();
      }
      
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  return new_plot;
}


bool ARTSGUI::MainMenu::Select(ArrayOfIndex& truths, const ArrayOfString& options, const String& menuname)
{
  bool pressed = false;
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu(menuname.c_str())) {
      for (Index i=0; i<truths.nelem(); i++) {
        bool val = truths[i];
        pressed |= ImGui::Checkbox(options[i].c_str(), &val);
        truths[i] = val;
        ImGui::Separator();
      }
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  return pressed;
}

bool ARTSGUI::MainMenu::Select(ArrayOfArrayOfIndex& truths, const ArrayOfString& dropdowns, const ArrayOfArrayOfString& options, const String& menuname)
{
  bool pressed = false;
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu(menuname.c_str())) {
      for (Index i=0; i<truths.nelem(); i++) {
        if (ImGui::BeginMenu(dropdowns[i].c_str())) {
          for (Index j=0; j<truths[i].nelem(); j++) {
            bool val = truths[i].at(j);
            pressed |= ImGui::Checkbox(options[i][j].c_str(), &val);
            truths[i].at(j) = val;
            ImGui::Separator();
          }
          ImGui::EndMenu();
        }
        ImGui::Separator();
      }
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  return pressed;
}
