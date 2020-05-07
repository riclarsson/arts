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


bool ARTSGUI::MainMenu::SelectAtmosphere(Numeric& p, Numeric& t, VectorView vmrs, const ArrayOfString& spec_list)
{
  bool new_plot=false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Atmospheric properties")) {
      float temperature = float(t);
      if (ImGui::SliderFloat("Temperature [K]", &temperature, 150.0f, 600.0f)) {
        if (temperature < 0.01)
          t = 0.01; 
        else
          t = Numeric(temperature);
        new_plot = true;
      }
      ImGui::Separator();
      
      float pressure = float(p);
      if (ImGui::SliderFloat("Pressure [Pa]", &pressure, 0.01f, 1'000'000.0f, "%g", 10.0f)) {
        if (pressure < 0.01)
          p = 0.01; 
        else
          p = Numeric(pressure);
        new_plot = true;
      }
      ImGui::Separator();
      
      for (Index i=0; i<spec_list.nelem(); i++) {
        float vmr = float(vmrs[i] * 1e6);
        if (ImGui::SliderFloat(spec_list[i].c_str(), &vmr, 0.00f, 1'000'000.0f, "%g", 10.0f)) {
          if (vmr < 0)
            vmrs[i] = 0;
          else
            vmrs[i] = Numeric(vmr) / 1e6;
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

bool ARTSGUI::MainMenu::SelectLOS(VectorView los)
{
  bool pressed = false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Line-of-sight")) {
      for (Index i=0; i<los.nelem(); i++) {
        if (i == 0) {
          float za = float(los[0]);
          if (ImGui::SliderFloat("Zenith Angle", &za, 0.0f, 180.0f, "%.3f deg")) {
            if (za<0.0f) za = 0.0f;
            else if (za>180.0f) za = 180.0f;
            los[0] = Numeric(za);
            pressed = true;
          }
        } else if (i == 1) {
          float aa = float(los[1]);
          if (ImGui::SliderFloat("Azimuth Angle", &aa, 0.0f, 360.0f, "%.3f deg")) {
            if (aa<0.0f) aa = 0.0f;
            else if (aa>360.0f) aa = 360.0f;
            los[1] = Numeric(aa);
            pressed = true;
          }
        }
        ImGui::Separator();
      }
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  return pressed;
}

bool ARTSGUI::MainMenu::SelectMAG(VectorView mag)
{
  bool pressed = false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Magnetic Field")) {
      for (Index i=0; i<mag.nelem(); i++) {
        if (i == 0) {
          float u = float(mag[0]*1e9);
          if (ImGui::SliderFloat("u", &u, -500'000.0f, 500'000.0f, "%.3f nT")) {
            mag[0] = 1e-9 * Numeric(u);
            pressed = true;
          }
        } else if (i == 1) {
          float v = float(mag[1]*1e9);
          if (ImGui::SliderFloat("v", &v, -500'000.0f, 500'000.0f, "%.3f nT")) {
            mag[1] = 1e-9 * Numeric(v);
            pressed = true;
          }
        } else if (i == 2) {
          float w = float(mag[2]*1e9);
          if (ImGui::SliderFloat("w", &w, -500'000.0f, 500'000.0f, "%.3f nT")) {
            mag[2] = 1e-9 * Numeric(w);
            pressed = true;
          }
        }
        ImGui::Separator();
      }
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  return pressed;
}
