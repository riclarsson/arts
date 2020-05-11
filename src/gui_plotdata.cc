#include "gui_inc.h"
#include "gui_plotdata.h"
#include "gui_plot.h"


bool ARTSGUI::Plotting::PlotFrame(Frame& frame, Config& cfg, bool select_frequency)
{
  bool new_plot=false;
  if (ImGui::BeginPlot(frame.title().c_str(), frame.xlabel().c_str(), frame.ylabel().c_str(), {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
    for (auto& line: frame)
      ImGui::Plot(line.name().c_str(), line.getter(), (void*)&line, line.size());
    
    if (select_frequency and frame.nelem())
      new_plot |= ARTSGUI::PlotMenu::SelectFrequency(frame[0].xview(), frame, cfg);
    ImGui::EndPlot();
  }
  return new_plot;
}
