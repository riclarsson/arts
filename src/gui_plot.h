#ifndef GUI_PLOT_H
#define GUI_PLOT_H

#include "gui_inc.h"

// Extra includes
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

#include "gui_help.h"
#include "gui_plotdata.h"
#include "gui_menubar.h"

/** All gui related things; should subnamespace things */
namespace ARTSGUI {

/** All plot-menu items go in here */
namespace PlotMenu {

/** Plot range modifier 
 *
 * Modifies the frame range copy as required.
 * 
 * @param[in,out] frame A plotting frame
 */
void range(ARTSGUI::Plotting::Frame& frame);

/** Plot scale modifier
 * 
 * Modifies the lines X-Y calculations as required
 * 
 * @param[in] frame A plotting frame
 */
void scale(ARTSGUI::Plotting::Frame& frame);

bool SelectFrequency(ARTSGUI::Plotting::Data& f, VectorView f_grid, ARTSGUI::Plotting::Frame& frame, ARTSGUI::Config& cfg);
};  // PlotMenuGUI
};  // ARTSGUI
  

#endif  // GUI_PLOT_H
