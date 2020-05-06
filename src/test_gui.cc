#include "gui_menubar.h"
#include "gui_plot.h"

int main(int, char**)
{
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::MainMenu::Config config;
  
  // Test data
  constexpr int n=3000;
  auto x = ARTSGUI::Plotting::Data(n);
  auto y1 = ARTSGUI::Plotting::Data(n);
  auto y2 = ARTSGUI::Plotting::Data(n);
  Vector tmp(n);
  for (int i=0; i<n; i++)
    tmp[i] = (6.28*5*i)/n;
  x.set(tmp);
  
  for (int i=0; i<n; i++)
    tmp[i] = std::sin(tmp[i]);
  y1.set(tmp);
  
  for (int i=0; i<n; i++)
    tmp[i] = (6.28*5*i)/n;
  for (int i=0; i<n; i++)
    tmp[i] = std::cos(tmp[i]);
  y2.set(tmp);
  
  ARTSGUI::Plotting::Line line1("Sine-curve", &x, &y1);
  ARTSGUI::Plotting::Line line2("Cosine-curve", &x, &y2);
  
  ARTSGUI::Plotting::Frame frame("Test plot", "X", "Y");
  frame.push_back(line1);
  frame.push_back(line2);
  
  // Style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
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
    if (ImGui::BeginPlot(frame.title().c_str(), frame.xlabel().c_str(), frame.ylabel().c_str(), {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
      for (auto& line: frame)
        ImGui::Plot(line.name().c_str(), line.getter(), (void*)&line, line.size());
      ImGui::EndPlot();
    }
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(frame);
    ARTSGUI::PlotMenu::scale(frame);
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  ARTSGUI::MainMenu::imgui_help(config);

  EndWhileLoopARTSGUI;
  CleanupARTSGUI;

  return 0;
}
