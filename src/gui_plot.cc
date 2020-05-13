#include "gui_plot.h"


void ARTSGUI::PlotMenu::range(ARTSGUI::Plotting::Frame& frame) 
{
  bool enable_paste = not frame.invalid_range();
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Frames")) {
      if (ImGui::BeginMenu(frame.title().c_str())) {
        ImPlotLimits limits = ImGui::GetPlotLimits();
        bool change = false;
        change |= ImGui::InputFloat2("X", &limits.X.Min, "%g", ImGuiInputTextFlags_CharsScientific);
        change |= ImGui::InputFloat2("Y", &limits.Y.Min, "%g", ImGuiInputTextFlags_CharsScientific);
        
        ImGui::Separator();
        
        if (ImGui::MenuItem("Copy"))
          frame.limits(limits);
        
        if (ImGui::MenuItem("Paste", NULL, false, enable_paste)) {
          limits = frame.limits();
          change = true;
        }
        
        if (change) ImGui::SetPlotLimits(limits);
        
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

void change_f(VectorView f, Numeric f0, Numeric f1, Numeric scale, Numeric offset)
{
  if (f0 < (1.0 - offset) / scale)
    f0 = (1.0 - offset) / scale;
  if (f1 < f0 + (1.0 - offset) / scale)
    f1 = f0 + (1.0 - offset) / scale;
  
  
  Numeric step = scale * (f1 - f0) / (Numeric(f.nelem()) - 1);
  for (Index i = 0; i < f.nelem() - 1; i++)
    f[i] = scale * f0 + Numeric(i) * step + offset;
  f[f.nelem() - 1] = scale * f1 + offset;
}

bool ARTSGUI::PlotMenu::SelectFrequency(ARTSGUI::Plotting::Frame& frame, ARTSGUI::Config& cfg)
{
  bool new_plot = false;
  
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Frequency")) {
      
      ImGui::Checkbox("Autoscale", &cfg.autoscale_x);
      ImGui::Separator();
      
      for (auto& line: frame) {
        if (ImGui::BeginMenu("Select Frequency Range")) {
          if (ImGui::MenuItem("Hz")) {line.x() -> scale(1.0); frame.xlabel("Frequency [Hz]");}
          if (ImGui::MenuItem("kHz")) {line.x() -> scale(1'000.0); frame.xlabel("Frequency [kHz]");}
          if (ImGui::MenuItem("MHz")) {line.x() -> scale(1'000'000.0); frame.xlabel("Frequency [MHz]");}
          if (ImGui::MenuItem("GHz")) {line.x() -> scale(1'000'000'000.0); frame.xlabel("Frequency [GHz]");}
          if (ImGui::MenuItem("THz")) {line.x() -> scale(1'000'000'000'000.0); frame.xlabel("Frequency [THz]");}
          if (ImGui::MenuItem("PHz")) {line.x() -> scale(1'000'000'000'000'000.0); frame.xlabel("Frequency [PHz]");}
          if (ImGui::MenuItem("EHz")) {line.x() -> scale(1'000'000'000'000'000'000.0); frame.xlabel("Frequency [EHz]");}
          ImGui::EndMenu();
        }
        ImGui::Separator();
        
        if (ImGui::BeginMenu((frame.title() + ": " + line.name()).c_str(), not cfg.autoscale_x))  {
          Numeric f0 = line.x() -> get(0);
          Numeric f1 = line.x() -> get(int(line.x()->view().nelem()-1));
          const Numeric negoffset = (1.0 - line.x()->offset()) / line.x()->scale();
          const Numeric posoffset = 10'000 - line.x()->offset() / line.x()->scale();
          if (ImGui::SliderScalar(String("Min " + frame.xlabel()).c_str(), ImGuiDataType_Double, &f0, &negoffset, &f1)) {
            change_f(line.x()->view(), f0, f1, line.x() -> scale(), line.x() -> offset());
            new_plot = true;
          }
          if (ImGui::SliderScalar(String("Max " + frame.xlabel()).c_str(), ImGuiDataType_Double, &f1, &f0, &posoffset)) {
            change_f(line.x()->view(), f0, f1, line.x() -> scale(), line.x() -> offset());
            new_plot = true;
          }
          ImGui::Separator();
          
          Numeric o = line.x() -> offset() / line.x() -> scale();
          const Numeric o0=-100.0, o1=100.0;
          if (ImGui::SliderScalar((String("Offset ") + frame.xlabel()).c_str(), ImGuiDataType_Double, &o, &o0, &o1)) {
            line.x() -> offset(o * line.x() -> scale());
          }
          ImGui::Separator();
          
          ImGui::EndMenu();
        }
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    
    if (cfg.autoscale_x) {
      if (frame.invalid_range()) {
        frame.limits(ImGui::GetPlotLimits());
      } else {
        ImPlotLimits limits = ImGui::GetPlotLimits();
        if (limits.X.Max not_eq frame.limits().X.Max or limits.X.Min not_eq frame.limits().X.Min) {
          frame.limits(limits);
          
          Numeric f0 = frame.limits().X.Min;
          Numeric f1 = frame.limits().X.Max;
          for (auto& line: frame)
            change_f(line.x() -> view(), f0, f1, line.x() -> scale(), line.x() -> offset());
          new_plot = true;
        }
      }
    }
    
    ImGui::EndMainMenuBar();
  }
  
  return new_plot;
}
