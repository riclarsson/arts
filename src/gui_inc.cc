#include "gui_inc.h"

void ARTSGUI::LayoutAndStyleSettings()
{
  auto& style = ImGui::GetStyle();
  style.FramePadding = {1.0f, 1.0f};
  style.WindowPadding = {1.0f, 1.0f};
  style.WindowRounding = 0.0f;
  style.WindowBorderSize = 0.0f;
}
