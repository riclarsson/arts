#include "gui_plot.h"
#include "gui_menubar.h"

int main(int, char**)
{
  InitializeARTSGUI;
  
  // Our states
  ImPlotRange copy_range;
  ARTSGUI::MainMenu::Config config;
  
  // Test data
  auto x = PlotData("X-axis", 50000);
  auto y = PlotData("Y-axis", 50000);
  Vector tmp(50000);
  for (int i=0; i<50000; i++)
    tmp[i] = (6.28*5*i)/100000;
  x.set(tmp);
  
  for (int i=0; i<50000; i++)
    tmp[i] = std::sin(tmp[i]);
  y.set(tmp);            
  
  Line line("XY-plot", &x, &y);
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  ARTSGUI::MainMenu::background(clear_color);
  ARTSGUI::MainMenu::imgui_help(config);
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  //Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 pos = ImGui::GetCursorPos();
  ImVec2 size = {float(width)-2*pos.x, float(height)-pos.y};
  
  // Show a simple window
  ImGui::SetNextWindowPos(pos);
  ImGui::SetNextWindowSize(size);
  if (ImGui::Begin("Plot tool", NULL, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
    if (ImGui::BeginPlot("XY-plot", "X-values", "Y-values", {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
      ImGui::Plot(line.name(), line.getter(), (void*)&line, line.size());
      ImGui::EndPlot();
    }
    
    // Menu bar
    ARTSGUI::PlotMenu::range(copy_range);
    ARTSGUI::PlotMenu::scale(line);
  }
  ImGui::End();

  EndWhileLoopARTSGUI;
  CleanupARTSGUI;

  return 0;
}
