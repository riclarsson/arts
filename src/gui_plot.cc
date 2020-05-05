#include "gui_plot.h"


void ARTSGUI::PlotMenu::range(ImPlotRange& copy) 
{
  bool enable_paste = not (std::isnan(copy.XMax) or std::isnan(copy.XMin) or 
  std::isnan(copy.YMax) or std::isnan(copy.YMin));
  
  if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("Range")) {
      ImPlotRange range = ImGui::GetPlotRange();
      bool change = false;
      change |= ImGui::InputFloat2("X", &range.XMin, "%g", ImGuiInputTextFlags_CharsScientific);
      change |= ImGui::InputFloat2("Y", &range.YMin, "%g", ImGuiInputTextFlags_CharsScientific);
      
      ImGui::Separator();
      
      if (ImGui::MenuItem("Copy"))
        copy = range;
      
      
      if (ImGui::MenuItem("Paste", NULL, false, enable_paste)) {
        range = copy;
        change = true;
      }
      
      if (change) ImGui::SetPlotRange(range);
      ImGui::EndMenu();
    }
  } 
  ImGui::EndMenuBar();
}

void ARTSGUI::PlotMenu::scale(Line& line)
{
  int N = int(line.TypeModifier());
  if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("Scale")) {
      
      
      if (ImGui::SliderInt("Linear Average", &N, 1, int(line.maxsize()))) {
        line.Linear(N);
      }
      
      
      ImGui::EndMenu();
    }
  }
  ImGui::EndMenuBar();
}
