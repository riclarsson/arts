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

bool ARTSGUI::PlotMenu::SelectFrequency(ARTSGUI::Plotting::Data& f, VectorView f_grid, ARTSGUI::Plotting::Frame& frame, ARTSGUI::Config& cfg)
{
  bool new_plot = false;
  
  if (ImGui::BeginMainMenuBar()) {
    bool change  = false;
    float f0 = float(f_grid[0]);
    float f1 = float(f_grid[f_grid.nelem()-1]);
    
    if (ImGui::BeginMenu("Frequency")) {
      
      ImGui::Checkbox("Autoscale", &cfg.autoscale_x);
      
      if (ImGui::SliderFloat("Min Freq", &f0, 0.01f, f1, "%g", 10.0f)) change = true;
      ImGui::Separator();
      if (ImGui::SliderFloat("Max Freq", &f1, f0, 1e15f, "%g", 10.0f)) change = true;
      ImGui::Separator();
      ImGui::EndMenu();
    }
    
    if (cfg.autoscale_x) {
      if (frame.invalid_range())
        frame.range(ImGui::GetPlotRange());
      else {
        ImPlotRange range = ImGui::GetPlotRange();
        if (range.XMax not_eq frame.range().XMax or range.XMin not_eq frame.range().XMin) {
          frame.range(range);
          change = true;
        }
      }
      f0 = frame.range().XMin;
      f1 = frame.range().XMax;
    }
    
    if (change) {
      if (f0 < 0.01f)
        f0 = 0.01f;
      if (f1 < f0 + 0.01f)
        f1 = f0 + 0.01f;
      
      Numeric step = Numeric(f1 - f0) / (Numeric(f_grid.nelem()) - 1);
      for (Index i = 0; i < f_grid.nelem() - 1; i++) f_grid[i] = Numeric(f0) + Numeric(i) * step;
      f_grid[f_grid.nelem() - 1] = Numeric(f1);
      f.set(f_grid);
      new_plot = true;
    }
    
    ImGui::EndMainMenuBar();
  }
  
  return new_plot;
}
