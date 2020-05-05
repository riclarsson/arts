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

/** All gui related things; should subnamespace things */
namespace ARTSGUI {

/** All plot-menu items go in here */
namespace PlotMenu {

/** Plot range modifier 
 *
 * @param[in] copy A copy of a range, if anyone is nan, the paste option is disabled
 */
void range(ImPlotRange& copy);

/** Plot scale modifier
 * 
 * @param[in] line A description of a line whose scales can be modified
 */
void scale(Line& line);
};  // PlotMenuGUI
};  // ARTSGUI
  

#endif  // GUI_PLOT_H
