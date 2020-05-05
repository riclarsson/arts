#include "gui_menubar.h"
#include "gui_help.h"

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

void ARTSGUI::MainMenu::background(ImVec4& colors)
{
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("View")) {
      ImGui::SliderFloat("Red", &colors.x, 0.0f, 1.0f);
      ImGui::SliderFloat("Green", &colors.y, 0.0f, 1.0f);
      ImGui::SliderFloat("Blue", &colors.z, 0.0f, 1.0f);
      ImGui::SliderFloat("Alpha", &colors.w, 0.0f, 1.0f);
      ImGui::Separator();
      ImGui::EndMenu();
    }
  ImGui::EndMainMenuBar();
  }
}
