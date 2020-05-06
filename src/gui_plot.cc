#include "gui_plot.h"


void ARTSGUI::PlotMenu::range(ARTSGUI::Plotting::Frame& frame) 
{
  bool enable_paste = not frame.invalid_range();
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Frames")) {
      if (ImGui::BeginMenu(frame.title().c_str())) {
        ImPlotRange range = ImGui::GetPlotRange();
        bool change = false;
        change |= ImGui::InputFloat2("X", &range.XMin, "%g", ImGuiInputTextFlags_CharsScientific);
        change |= ImGui::InputFloat2("Y", &range.YMin, "%g", ImGuiInputTextFlags_CharsScientific);
        
        ImGui::Separator();
        
        if (ImGui::MenuItem("Copy"))
          frame.range(range);
        
        if (ImGui::MenuItem("Paste", NULL, false, enable_paste)) {
          range = frame.range();
          change = true;
        }
        
        if (change) ImGui::SetPlotRange(range);
        
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
}

void ARTSGUI::PlotMenu::scale(ARTSGUI::Plotting::Frame& frame)
{
  for (auto& line: frame) {
    int N = int(line.TypeModifier());
    
    if (ImGui::BeginMainMenuBar()) {
      if (ImGui::BeginMenu("Frames")) {
        if (ImGui::BeginMenu(frame.title().c_str())) {
          if (ImGui::BeginMenu(line.name().c_str())) {
            if (ImGui::BeginMenu("Scale")) {
              if (ImGui::MenuItem("Reset Scale")) line.Linear(1);
              if (ImGui::SliderInt("Linear Average (quick)", &N, 1, int(line.maxsize()) / 2)) line.Linear(N);
              if (ImGui::SliderInt("Running Average (slow)", &N, 1, int(line.maxsize()) / 2)) line.RunningLinear(N);
              ImGui::Separator();
              ImGui::EndMenu();
            }
            ImGui::Separator();
            ImGui::EndMenu();
          }
          ImGui::Separator();
          ImGui::EndMenu();
        }
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::EndMainMenuBar();
    }
  }
}
