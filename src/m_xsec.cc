#include "absorption.h"
#include "auto_md.h"
#include "gui_menubar.h"
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
  new_plot |= ARTSGUI::MainMenu::SelectAtmosphere(abs_p, abs_t, abs_vmrs, vmr_list);
  
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

void PlotSpeciesLinesXsec(
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
  new_plot |= ARTSGUI::MainMenu::SelectAtmosphere(abs_p, abs_t, abs_vmrs, vmr_list);
  
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
        xsec = 0;
        src = 0;
        for (auto& dx: dxsec) dx = 0;
        for (auto& ds: dsrc) ds = 0;
        
        for (auto& lines: abs_lines_per_species[species]) {
          xsec_species(xsec, src, dummy1, dxsec, dsrc, dummy2,
            jacobian_quantities, jac_pos, f_grid, abs_p, abs_t,
            rtp_nlte, abs_vmrs, abs_species, lines,
            isotopologue_ratios.getIsotopologueRatio(lines.QuantumIdentity()),
                      partition_functions.getParamType(lines.QuantumIdentity()),
                      partition_functions.getParam(lines.QuantumIdentity()));
        }
        
        xsec_data[iline].overwrite(xsec(joker, 0));
        xsec_lines[iline] = ARTSGUI::Plotting::Line(species_list[species], &f, &xsec_data[iline]);
        iline += 1;
      }
    }
    
    xsec_frame.lines(xsec_lines);
  }
  
  EndWhileLoopARTSGUI;
  CleanupARTSGUI;
}  
