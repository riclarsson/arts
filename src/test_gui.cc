#include "gui_menubar.h"
#include "gui_plot.h"
#include "gui_windows.h"
#include "constants.h"

int main(int, char**)
{
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::Config config;
  
  // Test data
  constexpr int n=300;
  constexpr Numeric fac = 1.0 / Numeric(2*n);
  auto x = ARTSGUI::Plotting::Data(n);
  auto y1 = ARTSGUI::Plotting::Data(n);
  auto y2 = ARTSGUI::Plotting::Data(n);
  auto y3 = ARTSGUI::Plotting::Data(n);
  Vector tmp(n);
  for (int i=0; i<n; i++)
    tmp[i] =  -Constant::pi/2 + (Constant::pi*i)/n + fac;
  tmp[n-1] = Constant::pi/2 - fac;
  x.set(tmp);
  
  for (int i=0; i<n; i++)
    tmp[i] = std::sin(-Constant::pi/2 + (Constant::pi*i)/n + fac);
  tmp[n-1] = std::sin(Constant::pi/2 - fac);
  y1.set(tmp);
  for (int i=0; i<n; i++)
    tmp[i] = std::cos(-Constant::pi/2 + (Constant::pi*i)/n + fac);
  tmp[n-1] = std::cos(Constant::pi/2 - fac);
  y2.set(tmp);
  for (int i=0; i<n; i++)
    tmp[i] = std::tan(-Constant::pi/2 + (Constant::pi*i)/n + fac);
  tmp[n-1] = std::tan(Constant::pi/2 - fac);
  y3.set(tmp);
  
  ARTSGUI::Plotting::Line line1("Sine-curve", &x, &y1);
  ARTSGUI::Plotting::Line line2("Cosine-curve", &x, &y2);
  ARTSGUI::Plotting::Line line3("Tangent-curve", &x, &y3);
  
  ARTSGUI::Plotting::Frame frame1("Test sine", "X", "Y", line1);
  ARTSGUI::Plotting::Frame frame2("Test cos", "X", "Y", line2);
  ARTSGUI::Plotting::Frame frame3("Test tan", "X", "Y", line3);
  
  // Style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  const auto pos = ImGui::GetCursorPos();
  
  if (ARTSGUI::Windows::sub<2, 2, 0, 0>(window, pos, "Plot tool 1")) {
    ARTSGUI::Plotting::PlotFrame(frame1, config, false);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(frame1);
    ARTSGUI::PlotMenu::scale(frame1);
  }
  ImGui::End();
  
  if (ARTSGUI::Windows::sub<2, 2, 0, 1>(window, pos, "Plot tool 2")) {
    ARTSGUI::Plotting::PlotFrame(frame2, config, false);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(frame2);
    ARTSGUI::PlotMenu::scale(frame2);
  }
  ImGui::End();
  
  if (ARTSGUI::Windows::sub<2, 3, 1, 1, 1, 2>(window, pos, "Plot tool 3")) {
    ARTSGUI::Plotting::PlotFrame(frame3, config, false);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(frame3);
    ARTSGUI::PlotMenu::scale(frame3);
  }
  ImGui::End();
  
  if (ARTSGUI::Windows::sub<2, 3, 1, 0>(window, pos, "Text tool 4")) {
    ImGui::Text("The layout can be somewhat crazy but it must be set at compile time...\n"
                "Note that each call to a windowing function takes a name and if that\n"
                "name is not unique, there is a hidden error in the layout...\n");
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  ARTSGUI::MainMenu::imgui_help(config);

  EndWhileLoopARTSGUI;
  CleanupARTSGUI;

  return 0;
}
