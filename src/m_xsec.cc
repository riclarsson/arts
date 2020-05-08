#include "absorption.h"
#include "auto_md.h"
#include "gui_menubar.h"
#include "gui_windows.h"
#include "gui_plot.h"

void PlotBandXsec(
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
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
  const Numeric& fmin,
  const Numeric& fmax,
  const Index& fnum,
  const Verbosity& verbosity)
{
  // Output of methods
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx;
  ArrayOfArrayOfMatrix dsrc_xsec_per_species_dx;
  ArrayOfIndex asa(abs_species.nelem());
  std::iota(asa.begin(), asa.end(), 0);
  
  // Modifiable variables
  Vector abs_p;
  Vector abs_t;
  Matrix abs_vmrs;
  Vector f_grid;
  AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, rtp_pressure, rtp_temperature, rtp_vmr, verbosity);
  nlinspace(f_grid, fmin, fmax, fnum);
  
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
  Matrix xsec(abs_xsec_per_species[0]);
  Matrix src = do_lte ? dummy1 : src_xsec_per_species[0];
  ArrayOfMatrix dxsec = do_jac ? dabs_xsec_per_species_dx[0] : ArrayOfMatrix(0);
  ArrayOfMatrix dsrc = (do_jac and not do_lte) ? dsrc_xsec_per_species_dx[0] : ArrayOfMatrix(0);
  
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
  ArrayOfArrayOfIndex check_species;
  for (Index i=0; i<species_list.nelem(); i++) {
    band_list.push_back(ArrayOfString());
    
    for (Index j=0; j<abs_lines_per_species[i].nelem(); j++) {
      band_list.back().push_back(species_list[i] + " band " + std::to_string(j+1));
    }
    
    check_species.push_back(ArrayOfIndex(band_list.back().nelem(), 0));
  }
  ArrayOfString vmr_list;
  for (Index i=0; i<species_list.nelem(); i++) {
    vmr_list.push_back(species_list[i] + " VMR [ppmv]");
  }
  
  // Plotting data
  ARTSGUI::Plotting::Data f(f_grid);
  ARTSGUI::Plotting::ArrayOfData xsec_data(0);
  ARTSGUI::Plotting::ArrayOfLine xsec_lines(0);
  ARTSGUI::Plotting::Frame xsec_frame("Cross-section", "Frequency [Hz]", "Cross-section [1/m]", xsec_lines);
  
  if (not is_unique(species_list))
    throw std::runtime_error("Requires a unique set of species to work");
  
  // Create the plotting window
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::Config config;
  
  // Our style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  bool new_plot = ARTSGUI::MainMenu::Select(check_species, species_list, band_list, "Select Band");
  
  // Set atmospheric variables
  new_plot |= ARTSGUI::MainMenu::SelectAtmosphere(abs_p[0], abs_t[0], abs_vmrs(joker, 0), vmr_list);
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  // Draw a fullscreen plotting window
  if (ARTSGUI::Windows::full(window, ImGui::GetCursorPos(), "Plot tool")) {
    ARTSGUI::Plotting::PlotFrame(xsec_frame);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(xsec_frame);
    ARTSGUI::PlotMenu::scale(xsec_frame);
    new_plot |= ARTSGUI::PlotMenu::SelectFrequency(f, f_grid, xsec_frame, config);
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  
  if (new_plot) {
    Index numlines=0;
    for (Index species=0; species<band_list.nelem(); species++) {
      for (Index band=0; band<band_list[species].nelem(); band++) {
        if (check_species[species][band]) {
          numlines += 1;
        }
      }
    }
    
    xsec_data.resize(numlines);
    xsec_lines.resize(numlines);
    Index iline=0;
    
    for (Index species=0; species<band_list.nelem(); species++) {
      for (Index band=0; band<band_list[species].nelem(); band++) {
        if (check_species[species][band]) {
          xsec = 0;
          src = 0;
          for (auto& dx: dxsec) dx = 0;
          for (auto& ds: dsrc) ds = 0;
        
          auto& lines = abs_lines_per_species[species][band];
          xsec_species(xsec, src, dummy1, dxsec, dsrc, dummy2,
                       jacobian_quantities, jac_pos, f_grid, abs_p, abs_t,
                       rtp_nlte, abs_vmrs, abs_species, lines,
                       isotopologue_ratios.getIsotopologueRatio(lines.QuantumIdentity()),
                       partition_functions.getParamType(lines.QuantumIdentity()),
                       partition_functions.getParam(lines.QuantumIdentity()));
        
          xsec_data[iline].overwrite(xsec(joker, 0));
          xsec_lines[iline] = ARTSGUI::Plotting::Line(band_list[species][band], &f, &xsec_data[iline]);
          iline += 1;
        }
      }
    }
    
    xsec_frame.lines(xsec_lines);
  }
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}

void PlotXsecAgenda(
  Workspace& ws,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Numeric& rtp_pressure,
  const Numeric& rtp_temperature,
  const EnergyLevelMap& rtp_nlte,
  const Vector& rtp_vmr,
  const Agenda& abs_xsec_agenda,
  const Numeric& fmin,
  const Numeric& fmax,
  const Index& fnum,
  const Verbosity& verbosity)
{
  // Output of methods
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx;
  ArrayOfArrayOfMatrix dsrc_xsec_per_species_dx;
  ArrayOfIndex asa(abs_species.nelem());
  std::iota(asa.begin(), asa.end(), 0);
  
  // Modifiable variables
  Vector abs_p;
  Vector abs_t;
  Matrix abs_vmrs;
  Vector f_grid;
  AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, rtp_pressure, rtp_temperature, rtp_vmr, verbosity);
  nlinspace(f_grid, fmin, fmax, fnum);
  
  abs_xsec_agendaExecute(ws, abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                         abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_t, rtp_nlte, abs_vmrs, abs_xsec_agenda);
  
  if (abs_xsec_per_species.nelem() == 0)
    throw std::runtime_error("Cannot draw anything since xsec is zeroed");
  
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
  ArrayOfIndex check_species(species_list.nelem(), false);
  check_species[0] = true;
  
  // Plotting data
  ARTSGUI::Plotting::Data f(f_grid);
  ARTSGUI::Plotting::ArrayOfData xsec_data(1, ARTSGUI::Plotting::Data(abs_xsec_per_species[0](joker, 0)));
  ARTSGUI::Plotting::ArrayOfLine xsec_lines(1, ARTSGUI::Plotting::Line(species_list[0].c_str(), &f, &xsec_data[0]));
  ARTSGUI::Plotting::Frame xsec_frame("Cross-section", "Frequency [Hz]", "Cross-section [1/m]", xsec_lines);
  
  if (not is_unique(species_list))
    throw std::runtime_error("Requires a unique set of species to work");
  
  // Create the plotting window
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::Config config;
  
  // Our style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  bool new_plot = ARTSGUI::MainMenu::Select(check_species, species_list, "Select Species");
  
  // Set atmospheric variables
  new_plot |= ARTSGUI::MainMenu::SelectAtmosphere(abs_p[0], abs_t[0], abs_vmrs(joker, 0), vmr_list);
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  // Draw a fullscreen plotting window
  if (ARTSGUI::Windows::full(window, ImGui::GetCursorPos(), "Plot tool")) {
    ARTSGUI::Plotting::PlotFrame(xsec_frame);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(xsec_frame);
    ARTSGUI::PlotMenu::scale(xsec_frame);
    new_plot |= ARTSGUI::PlotMenu::SelectFrequency(f, f_grid, xsec_frame, config);
  }
  ImGui::End();
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  
  if (new_plot) {
    abs_xsec_agendaExecute(ws, abs_xsec_per_species, src_xsec_per_species, dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx,
                           abs_species, jacobian_quantities, asa, f_grid, abs_p, abs_t, rtp_nlte, abs_vmrs, abs_xsec_agenda);
    
    Index numlines=0;
    for (Index species=0; species<species_list.nelem(); species++) {
      if (check_species[species]) {
        numlines += 1;
      }
    }
    xsec_data.resize(numlines);
    xsec_lines.resize(numlines);
    Index iline=0;
    
    for (Index species=0; species<species_list.nelem(); species++) {
      if (check_species[species]) {
        xsec_data[iline].overwrite(abs_xsec_per_species[species](joker, 0));
        xsec_lines[iline] = ARTSGUI::Plotting::Line(species_list[species], &f, &xsec_data[iline]);
        iline += 1;
      }
    }
    
    xsec_frame.lines(xsec_lines);
  }
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}

void PlotPropmatAgenda(
  Workspace& ws,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Numeric& rtp_pressure_in,
  const Numeric& rtp_temperature_in,
  const EnergyLevelMap& rtp_nlte,
  const Vector& rtp_vmr_in,
  const Vector& rtp_mag_in,
  const Vector& rtp_los_in,
  const Agenda& propmat_clearsky_agenda,
  const Numeric& fmin,
  const Numeric& fmax,
  const Index& fnum,
  const Verbosity&)
{
  // Output of methods
  ArrayOfPropagationMatrix propmat_clearsky, dpropmat_clearsky_dx;
  ArrayOfStokesVector nlte_source, dnlte_dx_source, nlte_dx_dsource_dx;
  Vector f_grid;
  nlinspace(f_grid, fmin, fmax, fnum);
  Numeric rtp_pressure = rtp_pressure_in;
  Numeric rtp_temperature = rtp_temperature_in;
  Vector rtp_vmr = rtp_vmr_in;
  Vector rtp_mag = rtp_mag_in;
  Vector rtp_los = rtp_los_in;

  propmat_clearsky_agendaExecute(ws, propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dx_dsource_dx,
                                 jacobian_quantities, f_grid, rtp_mag, rtp_los, rtp_pressure, rtp_temperature, rtp_nlte,
                                 rtp_vmr, propmat_clearsky_agenda);
  
  if (propmat_clearsky.nelem() == 0)
    throw std::runtime_error("Cannot draw anything since propagation matrix is zeroed");
  
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
  ArrayOfIndex check_species(species_list.nelem(), false);
  check_species[0] = true;
  
  // Plotting data
  ARTSGUI::Plotting::Data f(f_grid);
  ARTSGUI::Plotting::ArrayOfData abs_data(1, ARTSGUI::Plotting::Data(propmat_clearsky[0].Kjj()));
  ARTSGUI::Plotting::ArrayOfLine abs_lines(1, ARTSGUI::Plotting::Line(species_list[0].c_str(), &f, &abs_data[0]));
  ARTSGUI::Plotting::Frame abs_frame("Absorption", "Frequency [Hz]", "Absorption [1/m]", abs_lines);
  
  if (not is_unique(species_list))
    throw std::runtime_error("Requires a unique set of species to work");
  
  Index component_option=0;
  const Index max_component=propmat_clearsky[0].NumberOfNeededVectors();
  
  // Create the plotting window
  InitializeARTSGUI;
  
  // Our states
  ARTSGUI::Config config;
  
  // Our style
  ARTSGUI::LayoutAndStyleSettings();
  
  // Main loop
  BeginWhileLoopARTSGUI;
  
  // Main menu bar
  ARTSGUI::MainMenu::fullscreen(config, window);
  ARTSGUI::MainMenu::quitscreen(config, window);
  
  bool new_plot = ARTSGUI::MainMenu::Select(check_species, species_list, "Select Species");
  
  // Set atmospheric variables
  new_plot |= ARTSGUI::MainMenu::SelectAtmosphere(rtp_pressure, rtp_temperature, rtp_vmr, vmr_list);
  new_plot |= ARTSGUI::MainMenu::SelectLOS(rtp_los);
  new_plot |= ARTSGUI::MainMenu::SelectMAG(rtp_mag, config);
  
  // Plotting defaults:
  ImGui::GetPlotStyle().LineWeight = 4;
  
  // Draw a fullscreen plotting window
  if (ARTSGUI::Windows::full(window, ImGui::GetCursorPos(), "Plot tool")) {
    ARTSGUI::Plotting::PlotFrame(abs_frame);
    
    // Menu bar for plot
    ARTSGUI::PlotMenu::range(abs_frame);
    ARTSGUI::PlotMenu::scale(abs_frame);
    new_plot |= ARTSGUI::PlotMenu::SelectFrequency(f, f_grid, abs_frame, config);
  }
  ImGui::End();
  
  // Special component menu bar
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Component")) {
      if (ImGui::MenuItem("I", NULL, false, 1<=max_component)) {component_option = 0; new_plot = true;}
      if (ImGui::MenuItem("Q", NULL, false, 2<=max_component)) {component_option = 1; new_plot = true;}
      if (ImGui::MenuItem("U", NULL, false, 3<=max_component)) {component_option = 2; new_plot = true;}
      if (ImGui::MenuItem("V", NULL, false, 5<=max_component)) {component_option = 4; new_plot = true;}
      if (ImGui::MenuItem("u", NULL, false, 4<=max_component)) {component_option = 3; new_plot = true;}
      if (ImGui::MenuItem("v", NULL, false, 6<=max_component)) {component_option = 5; new_plot = true;}
      if (ImGui::MenuItem("w", NULL, false, 7<=max_component)) {component_option = 6; new_plot = true;}
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (new_plot) {
    propmat_clearsky_agendaExecute(ws, propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dx_dsource_dx,
                                   jacobian_quantities, f_grid, rtp_mag, rtp_los, rtp_pressure, rtp_temperature, rtp_nlte,
                                   rtp_vmr, propmat_clearsky_agenda);;
    
    Index numlines=0;
    for (Index species=0; species<species_list.nelem(); species++) {
      if (check_species[species]) {
        numlines += 1;
      }
    }
    abs_data.resize(numlines);
    abs_lines.resize(numlines);
    Index iline=0;
    
    for (Index species=0; species<species_list.nelem(); species++) {
      if (check_species[species]) {
        if (component_option == 0) abs_data[iline].overwrite(propmat_clearsky[species].Kjj());
        if (component_option == 1) abs_data[iline].overwrite(propmat_clearsky[species].K12());
        if (component_option == 2) abs_data[iline].overwrite(propmat_clearsky[species].K13());
        if (component_option == 3) abs_data[iline].overwrite(propmat_clearsky[species].K23());
        if (component_option == 4) abs_data[iline].overwrite(propmat_clearsky[species].K14());
        if (component_option == 5) abs_data[iline].overwrite(propmat_clearsky[species].K24());
        if (component_option == 6) abs_data[iline].overwrite(propmat_clearsky[species].K34());
        abs_lines[iline] = ARTSGUI::Plotting::Line(species_list[species], &f, &abs_data[iline]);
        iline += 1;
      }
    }
    
    abs_frame.lines(abs_lines);
  }
  
  // Add help menu at end
  ARTSGUI::MainMenu::arts_help();
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}
