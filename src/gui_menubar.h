#ifndef GUI_MENUBAR_H
#define GUI_MENUBAR_H

#include "gui_inc.h"

#include "matpackII.h"
#include "mystring.h"
#include "config.h"

namespace ARTSGUI {
namespace MainMenu {

void fullscreen(Config& cfg, GLFWwindow* window);
void quitscreen(Config& cfg, GLFWwindow* window);
void imgui_help(Config& cfg);
void arts_help();
bool SelectAtmosphere(VectorView abs_p, VectorView abs_t, MatrixView abs_vmrs, const ArrayOfString& spec_list, const Index level=0);
bool Select(ArrayOfIndex& truths, const ArrayOfString& options, const String& menuname);
bool Select(ArrayOfArrayOfIndex& truths, const ArrayOfString& dropdowns, const ArrayOfArrayOfString& options, const String& menuname);

};  // MainMenu
};  // ARTSGUI

#endif  // GUI_MENUBAR_H
