#include "absorption.h"
#include "auto_md.h"
#include "gui_menubar.h"
#include "gui_plot.h"

void PlotBandXsec(
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Vector& f_grid,
  const Numeric& rtp_pressure,
  const Numeric& rtp_temperature,
  const EnergyLevelMap& rtp_nlte,
  const Vector& rtp_vmr,
  const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const SpeciesAuxData& isotopologue_ratios,
  const SpeciesAuxData& partition_functions,
  const Index& abs_xsec_agenda_checked,
  const Index& lbl_checked,
  const Index& nlte_do,
  const Verbosity& verbosity)
{
  // Output of methods
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx;
  ArrayOfArrayOfMatrix dsrc_xsec_per_species_dx;
  ArrayOfIndex asa(abs_species.nelem());
  std::iota(asa.begin(), asa.end(), 0);
  
  // Abs variables
  Vector abs_p;
  Vector abs_t;
  Matrix abs_vmrs;
  AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, rtp_pressure, rtp_temperature, rtp_vmr, verbosity);
  
  // Run once to have some output and to check the individual data
  abs_xsec_per_speciesInit(abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                           abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_xsec_agenda_checked,
                           nlte_do, verbosity);
  abs_xsec_per_speciesAddLines(abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                               abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_t, rtp_nlte, abs_vmrs,
                               abs_lines_per_species, isotopologue_ratios, partition_functions, lbl_checked, verbosity);
  
  if (abs_xsec_per_species.nelem() == 0)
    throw std::runtime_error("Cannot draw anything since xsec is zeroed");
  
  // Setup constant variables for the drawing loop
  const bool do_jac = supports_propmat_clearsky(jacobian_quantities);
  const bool do_lte = rtp_nlte.Data().empty();
  const ArrayOfIndex jac_pos = equivalent_propmattype_indexes(jacobian_quantities);
  
  // Uninteresting data
  static Matrix dummy1(0, 0);
  static ArrayOfMatrix dummy2(0);
  
  // Setup dynamic variables for the drawing loop
  Index band = 0;
  Index species = 0;
  Matrix xsec(abs_xsec_per_species[species]);
  Matrix src = do_lte ? dummy1 : src_xsec_per_species[species];
  ArrayOfMatrix dxsec = do_jac ? dabs_xsec_per_species_dx[species] : ArrayOfMatrix(0);
  ArrayOfMatrix dsrc = (do_jac and not do_lte) ? dsrc_xsec_per_species_dx[species] : ArrayOfMatrix(0);
  
  // Plotting data
  ARTSGUI::Plotting::Data f(f_grid);
  ARTSGUI::Plotting::Data xsec_data(abs_xsec_per_species[species](joker, 0));
  ARTSGUI::Plotting::Line xsec_line("Cross-section", &f, &xsec_data);
  ARTSGUI::Plotting::Frame xsec_frame("Cross-section", "Frequency [Hz]", "Cross-section [1/m]", xsec_line);
  
  // Menu data
  ArrayOfString species_list;
  for (auto& specs: abs_species) {
    species_list.push_back(String(""));
    for (Index ispec=0; ispec<specs.nelem(); ispec++) {
      species_list.back() += specs[ispec].Name();
      if (ispec+1 not_eq specs.nelem())
        species_list.back() += ',';
    }
  }
  ArrayOfArrayOfString band_list;
  for (Index i=0; i<species_list.nelem(); i++) {
    band_list.push_back(ArrayOfString());
    
    for (Index j=0; j<abs_lines_per_species[i].nelem(); j++) {
      band_list.back().push_back(species_list[i] + " band " + std::to_string(j+1) + " with " + std::to_string(abs_lines_per_species[i][j].NumLines()) + 
      " lines; " + std::to_string(abs_lines_per_species[i][j].AllLines().front().F0()) + " to " + std::to_string(abs_lines_per_species[i][j].AllLines().back().F0()) + " Hz");
    }
  }
  ArrayOfString vmr_list;
  for (Index i=0; i<species_list.nelem(); i++) {
    vmr_list.push_back(species_list[i] + " VMR [ppmv]");
  }
  
  if (not is_unique(species_list))
    throw std::runtime_error("Requires a unique set of species to work");
  
  // Create the plotting window
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::MainMenu::Config config;
  
  // Our style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Compute a new plot if something changes
  bool new_plot = false;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  // Select band
  if (ImGui::BeginMainMenuBar()) {
    for (Index i=0; i<species_list.nelem(); i++) {
      if (ImGui::BeginMenu("Select Band")) {
        if (ImGui::BeginMenu(species_list[i].c_str())) {
          for (Index j=0; j<band_list[i].nelem(); j++) {
            if (ImGui::MenuItem(band_list[i][j].c_str())) {
              species = i;
              band = j;
              new_plot = true;
            }
            ImGui::Separator();
          }
          ImGui::EndMenu();
        }
        ImGui::Separator();
        ImGui::EndMenu();
      }
    }
    ImGui::EndMainMenuBar();
  }
  
  // Set atmospheric variables
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Atmospheric properties")) {
      float temperature = float(abs_t[0]);
      if (ImGui::SliderFloat("Temperature [K]", &temperature, 150.0f, 600.0f)) {
        if (temperature < 0.01)
          abs_t[0] = 0.01; 
        else
          abs_t[0] = Numeric(temperature);
        new_plot = true;
      }
      ImGui::Separator();
      
      float pressure = float(abs_p[0]);
      if (ImGui::SliderFloat("Pressure [Pa]", &pressure, 0.01f, 1'000'000.0f, "%g", 10.0f)) {
        if (pressure < 0.01)
          abs_p[0] = 0.01; 
        else
          abs_p[0] = Numeric(pressure);
        new_plot = true;
      }
      ImGui::Separator();
      
      for (Index i=0; i<vmr_list.nelem(); i++) {
        float vmr = float(abs_vmrs(i, 0) * 1e6);
        if (ImGui::SliderFloat(vmr_list[i].c_str(), &vmr, 0.00f, 1'000'000.0f, "%g", 10.0f)) {
          if (vmr < 0)
            abs_vmrs(i, 0) = 0;
          else
            abs_vmrs(i, 0) = Numeric(vmr) / 1e6;
          new_plot = true;
        }
        ImGui::Separator();
      }
      
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  //Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 pos = ImGui::GetCursorPos();
  ImVec2 size = {float(width)-2*pos.x, float(height)-pos.y};
  
  // Show a simple window
  ImGui::SetNextWindowPos(pos);
  ImGui::SetNextWindowSize(size);
  if (ImGui::Begin("Plot tool", NULL, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
    if (ImGui::BeginPlot(xsec_frame.title().c_str(), xsec_frame.xlabel().c_str(), xsec_frame.ylabel().c_str(), {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
      for (auto& line: xsec_frame)
        ImGui::Plot(line.name().c_str(), line.getter(), (void*)&line, line.size());
      ImGui::EndPlot();
    }
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(xsec_frame);
    ARTSGUI::PlotMenu::scale(xsec_frame);
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  
  if (new_plot) {
    xsec = 0;
    src = 0;
    for (auto& dx: dxsec) dx = 0;
    for (auto& ds: dsrc) ds = 0;
    
    if (abs_lines_per_species.nelem() > species) {
      if (abs_lines_per_species[species].nelem() > band) {
        const AbsorptionLines& lines = abs_lines_per_species[species][band];
        xsec_species(
          xsec,
          src,
          dummy1,
          dxsec,
          dsrc,
          dummy2,
          jacobian_quantities,
          jac_pos,
          f_grid,
          abs_p,
          abs_t,
          rtp_nlte,
          abs_vmrs,
          abs_species,
          lines,
          isotopologue_ratios.getIsotopologueRatio(lines.QuantumIdentity()),
          partition_functions.getParamType(lines.QuantumIdentity()),
          partition_functions.getParam(lines.QuantumIdentity()));
        
        xsec_data.set(xsec(joker, 0));
        xsec_frame.pop_back();
        xsec_frame.push_back(xsec_line);
      }
    }
  }
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}

void PlotSpeciesLinesXsec(
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Vector& f_grid,
  const Numeric& rtp_pressure,
  const Numeric& rtp_temperature,
  const EnergyLevelMap& rtp_nlte,
  const Vector& rtp_vmr,
  const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const SpeciesAuxData& isotopologue_ratios,
  const SpeciesAuxData& partition_functions,
  const Index& abs_xsec_agenda_checked,
  const Index& lbl_checked,
  const Index& nlte_do,
  const Verbosity& verbosity)
{
  // Output of methods
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx;
  ArrayOfArrayOfMatrix dsrc_xsec_per_species_dx;
  ArrayOfIndex asa(abs_species.nelem());
  std::iota(asa.begin(), asa.end(), 0);
  
  // Abs variables
  Vector abs_p;
  Vector abs_t;
  Matrix abs_vmrs;
  AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, rtp_pressure, rtp_temperature, rtp_vmr, verbosity);
  
  // Run once to have some output and to check the individual data
  abs_xsec_per_speciesInit(abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                           abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_xsec_agenda_checked,
                           nlte_do, verbosity);
  abs_xsec_per_speciesAddLines(abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                               abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_t, rtp_nlte, abs_vmrs,
                               abs_lines_per_species, isotopologue_ratios, partition_functions, lbl_checked, verbosity);
  
  if (abs_xsec_per_species.nelem() == 0)
    throw std::runtime_error("Cannot draw anything since xsec is zeroed");
  
  // Setup constant variables for the drawing loop
  const bool do_jac = supports_propmat_clearsky(jacobian_quantities);
  const bool do_lte = rtp_nlte.Data().empty();
  const ArrayOfIndex jac_pos = equivalent_propmattype_indexes(jacobian_quantities);
  
  // Uninteresting data
  static Matrix dummy1(0, 0);
  static ArrayOfMatrix dummy2(0);
  
  // Setup dynamic variables for the drawing loop
  Index species = 0;
  Matrix xsec(abs_xsec_per_species[species]);
  Matrix src = do_lte ? dummy1 : src_xsec_per_species[species];
  ArrayOfMatrix dxsec = do_jac ? dabs_xsec_per_species_dx[species] : ArrayOfMatrix(0);
  ArrayOfMatrix dsrc = (do_jac and not do_lte) ? dsrc_xsec_per_species_dx[species] : ArrayOfMatrix(0);
  
  // Plotting data
  ARTSGUI::Plotting::Data f(f_grid);
  ARTSGUI::Plotting::Data xsec_data(abs_xsec_per_species[species](joker, 0));
  ARTSGUI::Plotting::Line xsec_line("Cross-section", &f, &xsec_data);
  ARTSGUI::Plotting::Frame xsec_frame("Cross-section", "Frequency [Hz]", "Cross-section [1/m]", xsec_line);
  
  // Menu data
  ArrayOfString species_list;
  for (auto& specs: abs_species) {
    species_list.push_back(String(""));
    for (Index ispec=0; ispec<specs.nelem(); ispec++) {
      species_list.back() += specs[ispec].Name();
      if (ispec+1 not_eq specs.nelem())
        species_list.back() += ',';
    }
  }
  ArrayOfString vmr_list;
  for (Index i=0; i<species_list.nelem(); i++) {
    vmr_list.push_back(species_list[i] + " VMR [ppmv]");
  }
  
  if (not is_unique(species_list))
    throw std::runtime_error("Requires a unique set of species to work");
  
  // Create the plotting window
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::MainMenu::Config config;
  
  // Our style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Compute a new plot if something changes
  bool new_plot = false;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  // Select band
  if (ImGui::BeginMainMenuBar()) {
    for (Index i=0; i<species_list.nelem(); i++) {
      if (ImGui::BeginMenu("Select Band")) {
        if (ImGui::MenuItem(species_list[i].c_str())) {
          species = i;
          new_plot = true;
        }
        ImGui::Separator();
        ImGui::EndMenu();
      }
    }
    ImGui::EndMainMenuBar();
  }
  
  // Set atmospheric variables
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Atmospheric properties")) {
      float temperature = float(abs_t[0]);
      if (ImGui::SliderFloat("Temperature [K]", &temperature, 150.0f, 600.0f)) {
        if (temperature < 0.01)
          abs_t[0] = 0.01; 
        else
          abs_t[0] = Numeric(temperature);
        new_plot = true;
      }
      ImGui::Separator();
      
      float pressure = float(abs_p[0]);
      if (ImGui::SliderFloat("Pressure [Pa]", &pressure, 0.01f, 1'000'000.0f, "%g", 10.0f)) {
        if (pressure < 0.01)
          abs_p[0] = 0.01; 
        else
          abs_p[0] = Numeric(pressure);
        new_plot = true;
      }
      ImGui::Separator();
      
      for (Index i=0; i<vmr_list.nelem(); i++) {
        float vmr = float(abs_vmrs(i, 0) * 1e6);
        if (ImGui::SliderFloat(vmr_list[i].c_str(), &vmr, 0.00f, 1'000'000.0f, "%g", 10.0f)) {
          if (vmr < 0)
            abs_vmrs(i, 0) = 0;
          else
            abs_vmrs(i, 0) = Numeric(vmr) / 1e6;
          new_plot = true;
        }
        ImGui::Separator();
      }
      
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  //Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 pos = ImGui::GetCursorPos();
  ImVec2 size = {float(width)-2*pos.x, float(height)-pos.y};
  
  // Show a simple window
  ImGui::SetNextWindowPos(pos);
  ImGui::SetNextWindowSize(size);
  if (ImGui::Begin("Plot tool", NULL, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
    if (ImGui::BeginPlot(xsec_frame.title().c_str(), xsec_frame.xlabel().c_str(), xsec_frame.ylabel().c_str(), {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
      for (auto& line: xsec_frame)
        ImGui::Plot(line.name().c_str(), line.getter(), (void*)&line, line.size());
      ImGui::EndPlot();
    }
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(xsec_frame);
    ARTSGUI::PlotMenu::scale(xsec_frame);
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  
  if (new_plot) {
    xsec = 0;
    src = 0;
    for (auto& dx: dxsec) dx = 0;
    for (auto& ds: dsrc) ds = 0;
    
    if (abs_lines_per_species.nelem() > species) {
      for (auto& lines: abs_lines_per_species[species]) {
        xsec_species(
          xsec,
          src,
          dummy1,
          dxsec,
          dsrc,
          dummy2,
          jacobian_quantities,
          jac_pos,
          f_grid,
          abs_p,
          abs_t,
          rtp_nlte,
          abs_vmrs,
          abs_species,
          lines,
          isotopologue_ratios.getIsotopologueRatio(lines.QuantumIdentity()),
                    partition_functions.getParamType(lines.QuantumIdentity()),
                    partition_functions.getParam(lines.QuantumIdentity()));
      }
       
      xsec_data.set(xsec(joker, 0));
      xsec_frame.pop_back();
      xsec_frame.push_back(xsec_line);
    }
  }
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}  
