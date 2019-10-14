/* Copyright (C) 2000-2012
   Stefan Buehler  <sbuehler@ltu.se>
   Axel von Engeln <engeln@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  \file   absorption.cc

  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Stefan Buehler and Axel von Engeln
*/

#include "absorption.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <map>
#include "arts.h"
#include "auto_md.h"
#include "file.h"
#include "interpolation_poly.h"
#include "linescaling.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"

#include "global_data.h"
#include "linefunctions.h"
#include "partial_derivatives.h"

/** Mapping of species auxiliary type names to SpeciesAuxData::AuxType enum */
static const char* SpeciesAuxTypeNames[] = {"NONE",
                                            "ISORATIO",  //Built-in type
                                            "ISOQUANTUM",
                                            "PART_TFIELD",
                                            "PART_COEFF",  //Built-in type
                                            "PART_COEFF_VIBROT"};

// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromCoeff_OldUnused(
    Numeric temperature) const {
  Numeric result = 0.;
  Numeric exponent = 1.;

  Vector::const_iterator it;

  //      cout << "T: " << temperature << "\n";
  for (it = mqcoeff.begin(); it != mqcoeff.end(); ++it) {
    result += *it * exponent;
    exponent *= temperature;
    //      cout << "it: " << it << "\n";
    //      cout << "res: " << result << ", exp: " << exponent << "\n";
  }
  return result;
}

// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromData_OldUnused(
    Numeric temperature) const {
  GridPosPoly gp;
  gridpos_poly(gp, mqcoeffgrid, temperature, mqcoeffinterporder);
  Vector itw;
  interpweights(itw, gp);
  return interp(itw, mqcoeff, gp);
}

void SpeciesAuxData::InitFromSpeciesData() {
  using global_data::species_data;

  mparams.resize(species_data.nelem());
  mparam_type.resize(species_data.nelem());

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    const Index niso = species_data[isp].Isotopologue().nelem();
    mparams[isp].resize(niso);
    mparam_type[isp].resize(niso);
    for (Index iso = 0; iso < niso; iso++) {
      mparams[isp][iso].resize(0);
      mparam_type[isp][iso] = SpeciesAuxData::AT_NONE;
    }
  }
}

void SpeciesAuxData::setParam(const Index species,
                              const Index isotopologue,
                              const AuxType auxtype,
                              const ArrayOfGriddedField1& auxdata) {
  mparam_type[species][isotopologue] = auxtype;
  mparams[species][isotopologue] = auxdata;
}

void SpeciesAuxData::setParam(const String& artstag,
                              const String& auxtype,
                              const ArrayOfGriddedField1& auxdata) {
  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if (!hinit) {
    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];
      for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
        SpecIsoMap indicies(i, j);
        String buf = sr.Name() + "-" + sr.Isotopologue()[j].Name();
        ArtsMap[buf] = indicies;
      }
    }
    hinit = true;
  }

  Index species;
  Index isotopologue;

  // ok, now for the cool index map:
  // is this arts identifier valid?
  const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artstag);
  if (i == ArtsMap.end()) {
    ostringstream os;
    os << "ARTS Tag: " << artstag << " is unknown.";
    throw runtime_error(os.str());
  }

  SpecIsoMap id = i->second;

  // Set mspecies:
  species = id.Speciesindex();

  // Set misotopologue:
  isotopologue = id.Isotopologueindex();

  Index this_auxtype = 0;

  while (this_auxtype < AT_FINAL_ENTRY &&
         auxtype != SpeciesAuxTypeNames[this_auxtype])
    this_auxtype++;

  if (this_auxtype != AT_FINAL_ENTRY) {
    setParam(species, isotopologue, (AuxType)this_auxtype, auxdata);
  } else {
    ostringstream os;
    os << "Unknown SpeciesAuxData type: " << auxtype;
    std::runtime_error(os.str());
  }
}

const ArrayOfGriddedField1& SpeciesAuxData::getParam(
    const Index species, const Index isotopologue) const {
  return mparams[species][isotopologue];
}

Numeric SpeciesAuxData::getIsotopologueRatio(const SpeciesTag& st) const {
  Numeric val = 1.0;
  if (st.Isotopologue() > -1)
    val = getParam(st.Species(), st.Isotopologue())[0].data[0];
  return val;
}

Numeric SpeciesAuxData::getIsotopologueRatio(const QuantumIdentifier& qid) const {
  return getParam(qid.Species(), qid.Isotopologue())[0].data[0];
}

String SpeciesAuxData::getTypeString(const Index species,
                                     const Index isotopologue) const {
  assert(mparam_type[species][isotopologue] < AT_FINAL_ENTRY);
  return SpeciesAuxTypeNames[mparam_type[species][isotopologue]];
}

bool SpeciesAuxData::ReadFromStream(String& artsid,
                                    istream& is,
                                    Index nparams,
                                    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if (!hinit) {
    out3 << "  ARTS index table:\n";

    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];

      for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
        SpecIsoMap indicies(i, j);
        String buf = sr.Name() + "-" + sr.Isotopologue()[j].Name();

        ArtsMap[buf] = indicies;

        // Print the generated data structures (for debugging):
        // The explicit conversion of Name to a c-String is
        // necessary, because setw does not work correctly for
        // stl Strings.
        const Index& i1 = ArtsMap[buf].Speciesindex();
        const Index& i2 = ArtsMap[buf].Isotopologueindex();

        out3 << "  Arts Identifier = " << buf << "   Species = " << setw(10)
             << setiosflags(ios::left) << species_data[i1].Name().c_str()
             << "iso = " << species_data[i1].Isotopologue()[i2].Name().c_str()
             << "\n";
      }
    }
    hinit = true;
  }

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return true;

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return true;

    // @ as first character marks catalogue entry
    char c;
    extract(c, line, 1);

    // check for empty line
    if (c == '@') {
      comment = false;
    }
  }

  // read the arts identifier String
  istringstream icecream(line);

  icecream >> artsid;

  if (artsid.length() != 0) {
    Index mspecies;
    Index misotopologue;

    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
    if (i == ArtsMap.end()) {
      ostringstream os;
      os << "ARTS Tag: " << artsid << " is unknown.";
      throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;

    // Set mspecies:
    mspecies = id.Speciesindex();

    // Set misotopologue:
    misotopologue = id.Isotopologueindex();

    ArrayOfGriddedField1 ratios;
    ratios.resize(1);
    // Extract accuracies:
    try {
      Numeric p = NAN;
      std::vector<Numeric> aux;
      for (Index ip = 0; ip < nparams; ip++) {
        icecream >> double_imanip() >> p;
        aux.push_back(p);
      }

      Vector grid;
      if (aux.size() > 1)
        nlinspace(grid, 1, (Numeric)aux.size(), aux.size());
      else
        grid = Vector(1, .1);

      ratios[0].set_grid(0, grid);
      ratios[0].data = aux;
      mparams[mspecies][misotopologue] = ratios;
    } catch (const runtime_error&) {
      throw runtime_error("Error reading SpeciesAuxData.");
    }
  }

  // That's it!
  return false;
}

void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& isoratios) {
  using global_data::species_data;

  // Check total number of species:
  if (species_data.nelem() != isoratios.nspecies()) {
    ostringstream os;
    os << "Number of species in SpeciesAuxData (" << isoratios.nspecies()
       << ") "
       << "does not fit builtin species data (" << species_data.nelem() << ").";
    throw runtime_error(os.str());
  }

  // For the selected species, we check all isotopes by looping over the
  // species data. (Trying to check only the isotopes actually used gets
  // quite complicated, actually, so we do the simple thing here.)

  // Loop over the absorption species:
  for (Index i = 0; i < abs_species.nelem(); i++) {
    // sp is the index of this species in the internal lookup table
    const Index sp = abs_species[i][0].Species();

    // Get handle on species data for this species:
    const SpeciesRecord& this_sd = species_data[sp];

    // Check number of isotopologues:
    if (this_sd.Isotopologue().nelem() != isoratios.nisotopologues(sp)) {
      ostringstream os;
      os << "Incorrect number of isotopologues in isotopologue data.\n"
         << "Species: " << this_sd.Name() << ".\n"
         << "Number of isotopes in SpeciesAuxData ("
         << isoratios.nisotopologues(sp) << ") "
         << "does not fit builtin species data ("
         << this_sd.Isotopologue().nelem() << ").";
      throw runtime_error(os.str());
    }

    for (Index iso = 0; iso < this_sd.Isotopologue().nelem(); ++iso) {
      // For "real" species (not representing continau) the isotopologue
      // ratio must not be NAN or below zero.
      if (!this_sd.Isotopologue()[iso].isContinuum()) {
        if (std::isnan(isoratios.getParam(sp, iso)[0].data[0]) ||
            isoratios.getParam(sp, iso)[0].data[0] < 0.) {
          ostringstream os;
          os << "Invalid isotopologue ratio.\n"
             << "Species: " << this_sd.Name() << "-"
             << this_sd.Isotopologue()[iso].Name() << "\n"
             << "Ratio:   " << isoratios.getParam(sp, iso)[0].data[0];
          throw runtime_error(os.str());
        }
      }
    }
  }
}

void checkPartitionFunctions(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& partfun) {
  using global_data::species_data;

  // Check total number of species:
  if (species_data.nelem() != partfun.nspecies()) {
    ostringstream os;
    os << "Number of species in SpeciesAuxData (" << partfun.nspecies() << ") "
       << "does not fit builtin species data (" << species_data.nelem() << ").";
    throw runtime_error(os.str());
  }

  // For the selected species, we check all isotopes by looping over the
  // species data. (Trying to check only the isotopes actually used gets
  // quite complicated, actually, so we do the simple thing here.)

  // Loop over the absorption species:
  for (Index i = 0; i < abs_species.nelem(); i++) {
    // sp is the index of this species in the internal lookup table
    const Index sp = abs_species[i][0].Species();

    // Get handle on species data for this species:
    const SpeciesRecord& this_sd = species_data[sp];

    // Check number of isotopologues:
    if (this_sd.Isotopologue().nelem() != partfun.nisotopologues(sp)) {
      ostringstream os;
      os << "Incorrect number of isotopologues in partition function data.\n"
         << "Species: " << this_sd.Name() << ".\n"
         << "Number of isotopes in SpeciesAuxData ("
         << partfun.nisotopologues(sp) << ") "
         << "does not fit builtin species data ("
         << this_sd.Isotopologue().nelem() << ").";
      throw runtime_error(os.str());
    }
  }
}

void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(
    SpeciesAuxData& sad) {
  using global_data::species_data;

  sad.InitFromSpeciesData();

  Vector grid(1, 1.);
  ArrayOfGriddedField1 ratios;
  ratios.resize(1);
  ratios[0].set_name("IsoRatios");
  ratios[0].set_grid_name(0, "Index");
  ratios[0].set_grid(0, grid);
  ratios[0].resize(1);

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem();
         iiso++) {
      ratios[0].data[0] = species_data[isp].Isotopologue()[iiso].Abundance();
      sad.setParam(isp, iiso, SpeciesAuxData::AT_ISOTOPOLOGUE_RATIO, ratios);
    }
  }
}

void fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(
    SpeciesAuxData& sad) {
  using global_data::species_data;

  sad.InitFromSpeciesData();

  ArrayOfGriddedField1 pfuncs;
  pfuncs.resize(2);
  pfuncs[0].set_name("PartitionFunction");
  pfuncs[0].set_grid_name(0, "Coeff");
  pfuncs[1].set_grid_name(0, "Temperature");

  ArrayOfString tgrid;
  tgrid.resize(2);
  tgrid[0] = "Tlower";
  tgrid[1] = "Tupper";

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem();
         iiso++) {
      Vector grid;
      const Vector& coeffs = species_data[isp].Isotopologue()[iiso].GetCoeff();

      assert(coeffs.nelem() >= 2);

      nlinspace(grid, 0, (Numeric)coeffs.nelem() - 1., coeffs.nelem());
      pfuncs[0].set_grid(0, grid);
      pfuncs[0].data = coeffs;

      const Vector& temp_range =
          species_data[isp].Isotopologue()[iiso].GetCoeffGrid();

      // Temperature data should either contain two Ts, lower and upper value of
      // the valid range or be empty
      assert(temp_range.nelem() == 0 || temp_range.nelem() == 2);

      if (temp_range.nelem() == 2) {
        pfuncs[1].set_grid(0, tgrid);
        pfuncs[1].data = temp_range;
      } else {
        pfuncs[1].resize(0);
        pfuncs[1].set_grid(0, ArrayOfString());
      }

      sad.setParam(
          isp, iiso, SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, pfuncs);
    }
  }
}

ostream& operator<<(ostream& os, const SpeciesRecord& sr) {
  for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
    os << sr.Name() << "-" << sr.Isotopologue()[j].Name() << "\n";
  }
  return os;
}

ostream& operator<<(ostream& os, const SpeciesAuxData& sad) {
  using global_data::species_data;
  for (Index sp = 0; sp < sad.nspecies(); sp++) {
    for (Index iso = 0; iso < sad.nisotopologues(sp); iso++) {
      os << species_name_from_species_index(sp) << "-"
         << global_data::species_data[sp].Isotopologue()[iso].Name();
      os << " " << sad.getTypeString(sp, iso) << std::endl;
      for (Index ip = 0; ip < sad.getParam(sp, iso).nelem(); ip++)
        os << "AuxData " << ip << " " << sad.getParam(sp, iso) << std::endl;
    }
  }

  return os;
}

/** Calculate line absorption cross sections for one tag group. All
    lines in the line list must belong to the same species. This must
    be ensured by abs_lines_per_speciesCreateFromLines, so it is only verified
    with assert. Also, the input vectors abs_p, and abs_t must all
    have the same dimension.

    This is mainly a copy of abs_species which is removed now, with
    the difference that the vmrs are removed from the absorption
    coefficient calculation. (the vmr is still used for the self
    broadening)

    Continua are not handled by this function, you have to call
    xsec_continuum_tag for those.

    \param[in,out] xsec_attenuation   Cross section of one tag group. This is now the
                                      true absorption cross section in units of m^2.
    \param[in,out] xsec_source        Cross section of one tag group. This is now the
                                      true source cross section in units of m^2.
    \param[in,out] xsec_phase         Cross section of one tag group. This is now the
                                      true dispersion cross section in units of m^2.
    \param[in] f_grid       Frequency grid.
    \param[in] abs_p        Pressure grid.
    \param[in] abs_t        Temperatures associated with abs_p.
    \param[in] all_vmrs     Gas volume mixing ratios [nspecies, np].
    \param[in] abs_species  Species tags for all species.
    \param[in] this_species Index of the current species in abs_species.
    \param[in] abs_lines    The spectroscopic line list.
    \param[in] ind_ls       Index to lineshape function.
    \param[in] ind_lsn      Index to lineshape norm.
    \param[in] cutoff       Lineshape cutoff.
    \param[in] isotopologue_ratios  Isotopologue ratios.

    \author Stefan Buehler and Axel von Engeln
    \date   2001-01-11 

    Changed from pseudo cross sections to true cross sections

    \author Stefan Buehler 
    \date   2007-08-08 

    Adapted to new Perrin line parameters, treating broadening by different 
    gases explicitly
 
    \author Stefan Buehler
    \date   2012-09-03
    
    Adapted to LineRecord independency in central parts of function.
    
    \author Richard Larsson
    \date   2014-10-29

*/
void xsec_species(MatrixView xsec_attenuation,
                  MatrixView xsec_source,
                  MatrixView xsec_phase,
                  ConstVectorView f_grid,
                  ConstVectorView abs_p,
                  ConstVectorView abs_t,
                  ConstMatrixView abs_t_nlte,
                  ConstMatrixView all_vmrs,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const ArrayOfLineRecord& abs_lines,
                  const Index ind_ls,
                  const Index ind_lsn,
                  const Numeric cutoff,
                  const SpeciesAuxData& isotopologue_ratios,
                  const SpeciesAuxData& partition_functions,
                  const Verbosity& verbosity) {
  CREATE_OUT3;

  extern const Numeric DOPPLER_CONST;
  static const Numeric doppler_const = DOPPLER_CONST;

  // dimension of f_grid, abs_lines
  const Index nf = f_grid.nelem();
  const Index nl = abs_lines.nelem();

  // number of pressure levels:
  const Index np = abs_p.nelem();

  // Define the vector for the line shape function and the
  // normalization factor of the lineshape here, so that we don't need
  // so many free store allocations.  the last element is used to
  // calculate the value at the cutoff frequency
  Vector ls_attenuation(nf + 1);
  Vector ls_phase_dummy;
  Vector fac(nf + 1);

  const bool cut = (cutoff != -1) ? true : false;

  const bool calc_phase = 0;

  // Check that the frequency grid is sorted in the case of lineshape
  // with cutoff. Duplicate frequency values are allowed.
  if (cut) {
    if (!is_sorted(f_grid)) {
      std::ostringstream os;
      os << "If you use a lineshape function with cutoff, your\n"
         << "frequency grid *f_grid* must be sorted.\n"
         << "(Duplicate values are allowed.)";
      throw std::runtime_error(os.str());
    }
  }

  // Check that all temperatures are non-negative
  bool negative = false;

  for (Index i = 0; !negative && i < abs_t.nelem(); i++) {
    if (abs_t[i] < 0.) negative = true;
  }

  if (negative) {
    std::ostringstream os;
    os << "abs_t contains at least one negative temperature value.\n"
       << "This is not allowed.";
    throw std::runtime_error(os.str());
  }

  // We need a local copy of f_grid which is 1 element longer, because
  // we append a possible cutoff to it.
  // The initialization of this has to be inside the line loop!
  Vector f_local(nf + 1);
  f_local[Range(0, nf)] = f_grid;

  // Voigt generally needs a different frequency grid. If we allocate
  // that in the outer loop, instead of in voigt, we don't have the
  // free store allocation at each lineshape call. Calculation is
  // still done in the voigt routine itself, this is just an auxillary
  // parameter, passed to lineshape. For selected lineshapes (e.g.,
  // Rosenkranz) it is used additionally to pass parameters needed in
  // the lineshape (e.g., overlap, ...). Consequently we have to
  // assure that aux has a dimension not less then the number of
  // parameters passed.
  Index ii = (nf + 1 < 10) ? 10 : nf + 1;
  Vector aux(ii);

  // Check that abs_p, abs_t, and abs_vmrs have consistent
  // dimensions. This could be a user error, so we throw a
  // runtime_error.

  if (abs_t.nelem() != np) {
    std::ostringstream os;
    os << "Variable abs_t must have the same dimension as abs_p.\n"
       << "abs_t.nelem() = " << abs_t.nelem() << '\n'
       << "abs_p.nelem() = " << np;
    throw std::runtime_error(os.str());
  }

  // all_vmrs should have dimensions [nspecies, np]:

  if (all_vmrs.ncols() != np) {
    std::ostringstream os;
    os << "Number of columns of all_vmrs must match abs_p.\n"
       << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
       << "abs_p.nelem() = " << np;
    throw std::runtime_error(os.str());
  }

  const Index nspecies = abs_species.nelem();

  if (all_vmrs.nrows() != nspecies) {
    std::ostringstream os;
    os << "Number of rows of all_vmrs must match abs_species.\n"
       << "all_vmrs.nrows() = " << all_vmrs.nrows() << '\n'
       << "abs_species.nelem() = " << nspecies;
    throw std::runtime_error(os.str());
  }

  // Check that the dimension of xsec is indeed [f_grid.nelem(),
  // abs_p.nelem()]:
  if (xsec_attenuation.nrows() != nf || xsec_attenuation.ncols() != np) {
    std::ostringstream os;
    os << "Variable xsec_attenuation must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
       << "[xsec_attenuation.nrows(),xsec_attenuation.ncols()] = ["
       << xsec_attenuation.nrows() << ", " << xsec_attenuation.ncols() << "]\n"
       << "f_grid.nelem() = " << nf << '\n'
       << "abs_p.nelem() = " << np;
    throw std::runtime_error(os.str());
  }

  // Check that the dimension of xsec source is [f_grid.nelem(), abs_p.nelem()]:
  bool calc_src;
  if (xsec_source.nrows() != nf || xsec_source.ncols() != np) {
    if (xsec_source.nrows() != 0 || xsec_source.ncols() != 0) {
      std::ostringstream os;
      os << "Variable xsec_source must have dimensions [f_grid.nelem(),abs_p.nelem()] or [0,0].\n"
         << "[xsec_source.nrows(),xsec_source.ncols()] = ["
         << xsec_source.nrows() << ", " << xsec_source.ncols() << "]\n"
         << "f_grid.nelem() = " << nf << '\n'
         << "abs_p.nelem() = " << np;
      throw std::runtime_error(os.str());
    } else {
      calc_src = false;
    }
  } else {
    calc_src = true;
  }

  if (xsec_phase.nrows() != 0 || xsec_phase.ncols() != 0) {
    std::ostringstream os;
    os << "Variable xsec_phase must have dimensions [0,0] here.\n"
       << "[xsec_phase.nrows(),xsec_phase.ncols()] = [" << xsec_phase.nrows()
       << ", " << xsec_phase.ncols() << "]\n"
       << "f_grid.nelem() = " << nf << '\n'
       << "abs_p.nelem() = " << np;
    throw std::runtime_error(os.str());
  }

  String fail_msg;
  bool failed = false;

  // Loop all pressures:
  if (np)
#pragma omp parallel for if (!arts_omp_in_parallel() &&        \
                             np >= arts_omp_get_max_threads()) \
    firstprivate(ls_attenuation, fac, f_local, aux)
    for (Index i = 0; i < np; ++i) {
      if (failed) continue;

      // holder when things can be empty
      Vector empty_vector(0);

      // Store input profile variables, this is perhaps slightly faster.
      const Numeric p_i = abs_p[i];
      const Numeric t_i = abs_t[i];
      ConstVectorView t_nlte_i = calc_src ? abs_t_nlte(joker, i) : empty_vector;

      //out3 << "  p = " << p_i << " Pa\n";

      // Calculate total number density from pressure and temperature.
      // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
      //      const Numeric n = p_i / BOLTZMAN_CONST / t_i;
      // This is not needed anymore, since we now calculate true cross
      // sections, which do not contain the n.

      // Get handle on xsec for this pressure level i.
      // Watch out! This is output, we have to be careful not to
      // introduce race conditions when writing to it.
      VectorView xsec_i_attenuation = xsec_attenuation(Range(joker), i);
      VectorView xsec_i_source =
          calc_src ? xsec_source(Range(joker), i) : empty_vector;

      //       if (omp_in_parallel())
      //         cout << "omp_in_parallel: true\n";
      //       else
      //         cout << "omp_in_parallel: false\n";

      // Prepare a variable that can be used by the individual LBL
      // threads to add up absorption:
      Index n_lbl_threads;
      if (arts_omp_in_parallel()) {
        // If we already are running parallel, then the LBL loop
        // will not be parallelized.
        n_lbl_threads = 1;
      } else {
        n_lbl_threads = arts_omp_get_max_threads();
      }
      Matrix xsec_accum_attenuation(
          n_lbl_threads, xsec_i_attenuation.nelem(), 0);
      Matrix xsec_accum_source(n_lbl_threads, xsec_i_source.nelem());
      if (calc_src) xsec_accum_source = 0.0;

      // Simple caching of partition function to avoid recalculating things.
      Numeric qt_cache = -1, qref_cache = -1;
      Index iso_cache = -1;
      Numeric line_t_cache = -1;

      ConstVectorView vmrs = all_vmrs(joker, i);

      // Loop all lines:
      if (nl) {
#pragma omp parallel for if (!arts_omp_in_parallel() &&        \
                             nl >= arts_omp_get_max_threads()) \
    firstprivate(ls_attenuation,                               \
                 fac,                                          \
                 f_local,                                      \
                 aux,                                          \
                 qt_cache,                                     \
                 qref_cache,                                   \
                 iso_cache,                                    \
                 line_t_cache)
        for (Index l = 0; l < nl; ++l) {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          //           if (omp_in_parallel())
          //             cout << "LBL: omp_in_parallel: true\n";
          //           else
          //             cout << "LBL: omp_in_parallel: false\n";

          // The try block here is necessary to correctly handle
          // exceptions inside the parallel region.
          try {
            const LineRecord& l_l = abs_lines[l];

            // Prepare pressure broadening parameters
            Numeric gamma_0, gamma_2, eta, df_0, df_2, f_VC;
            l_l.SetPressureBroadeningParameters(gamma_0,
                                                gamma_2,
                                                eta,
                                                df_0,
                                                df_2,
                                                f_VC,
                                                t_i,
                                                p_i,
                                                vmrs,
                                                abs_species);

            // Check the chache is the tempearture of the line and the isotope is the same to avoid recalculating the partition sum
            if (iso_cache != l_l.Isotopologue() || line_t_cache != l_l.Ti0()) {
              iso_cache = l_l.Isotopologue();
              line_t_cache = l_l.Ti0();
              qref_cache =
                  -1;  // no need to reset qt since it is done internally.
            }

            // Prepare line strength scaling
            Numeric partition_ratio, K1, K2, abs_nlte_ratio, src_nlte_ratio;
            GetLineScalingData(
                qt_cache,
                qref_cache,
                partition_ratio,
                K1,
                K2,
                abs_nlte_ratio,
                src_nlte_ratio,
                partition_functions.getParamType(l_l.Species(),
                                                 l_l.Isotopologue()),
                partition_functions.getParam(l_l.Species(), l_l.Isotopologue()),
                t_i,
                l_l.Ti0(),
                l_l.F(),
                l_l.Elow(),
                calc_src,
                l_l.Evlow(),
                l_l.Evupp(),
                l_l.NLTELowerIndex(),
                l_l.NLTEUpperIndex(),
                t_nlte_i);

            // Dopple broadening
            const Numeric sigma = l_l.F() * doppler_const *
                                  sqrt(t_i / l_l.IsotopologueData().Mass());

            const bool calc_partials = false;
            Vector da_dF, dp_dF, da_dP, dp_dP;

            Range this_f_range(0, 0);

            // Calculate line cross section
            xsec_single_line(  // OUTPUT
                xsec_accum_attenuation(arts_omp_get_thread_num(), joker),
                xsec_accum_source(arts_omp_get_thread_num(), joker),
                empty_vector,
                // HELPER:
                ls_attenuation,
                ls_phase_dummy,
                da_dF,
                dp_dF,
                da_dP,
                dp_dP,
                this_f_range,
                fac,
                aux,
                // FREQUENCY
                f_local,
                f_grid,
                nf,
                cutoff,
                l_l.F(),
                // LINE STRENGTH
                l_l.I0(),
                partition_ratio,
                K1 * K2,
                abs_nlte_ratio,
                src_nlte_ratio,
                isotopologue_ratios
                    .getParam(l_l.Species(), l_l.Isotopologue())[0]
                    .data[0],
                // ATMOSPHERIC TEMPERATURE
                t_i,
                // LINE SHAPE
                ind_ls,
                ind_lsn,
                // LINE BROADENING
                gamma_0,
                gamma_2,
                eta,
                df_0,
                df_2,
                sigma,
                f_VC,
                // LINE MIXING
                0,
                0,
                0,
                // FEATURE FLAGS
                cut,
                calc_phase,
                calc_partials,
                calc_src);

          }  // end of try block
          catch (const std::runtime_error& e) {
#pragma omp critical(xsec_species_fail)
            {
              fail_msg = e.what();
              failed = true;
            }
          }

        }  // end of parallel LBL loop
      }

      // Bail out if an error occurred in the LBL loop
      if (failed) continue;

      // Now we just have to add up all the rows of xsec_accum:
      for (Index j = 0; j < xsec_accum_attenuation.nrows(); ++j) {
        xsec_i_attenuation += xsec_accum_attenuation(j, Range(joker));
        if (calc_src) {
          xsec_i_source += xsec_accum_source(j, Range(joker));
        }
      }
    }  // end of parallel pressure loop

  if (failed)
    throw std::runtime_error("Run-time error in function: xsec_species\n" +
                             fail_msg);
}

/** 

   Calculate line absorption cross sections for one line at one layer
   and accumulates to the total phase and attenuation change.
   
   No dependency on LineRecord to increase speed for wrapper applications.
   
   \param[out] xsec_accum_attenuation   Cross section of one tag group. This is now the
                                    true absorption cross section in units of m^2.
                                    It has inputs of all previously calculated lines.
   \param[out] xsec_accum_source        Deviation of line cross section from LTE cross-section
   \param[out] xsec_accum_phase         Cross section of one tag group. This is now the
                                    true dispersion cross section in units of m^2.
                                    It has inputs of all previously calculated lines.
   \param[out] attenuation              Input only to increase speed.  Holds attenuation internally.
   \param[out] phase                    Input only to increase speed.  Holds phase internally.
   \param[out] fac                      Input only to increase speed.  Holds lineshape factor internally.
   \param[out] aux                      Input only to increase speed.  Holds f_grid factor internally.
   \param[out] f_local                  Input only to increase speed.  Holds f_grid internally.
   \param[in] f_grid,                   Frequency grid
   \param[in] nf,                       Number of frequencies to calculate
   \param[in] cutoff,                   Lineshape cutoff.
   \param[in] F0,                       Line center
   \param[in] intensity,                Line intensity
   \param[in] part_fct_ratio            Ratio of partition sums to atmospheric temperature
   \param[in] boltzmann_ratio           Ratio of Boltzmann statistics to atmospheric temperature
   \param[in] abs_nlte_ratio            Ratio of absorption intensity to LTE
   \param[in] src_nlte_ratio            Ratio of emission intesity to LTE
   \param[in] Isotopologue_Ratio,       Ratio of the isotopologue in the atmosphere
   \param[in] temperature,              Atmospheric temperature
   \param[in] ind_ls,                   Index to lineshape function.
   \param[in] ind_lsn,                  Index to lineshape norm.
   \param[in] gamma_0,                  Line pressure broadening
   \param[in] gamma_2,                  Speed-dependent line pressure broadening
   \param[in] eta,                      Correlation of pressure broadening parameters
   \param[in] df_0,                     Line pressure shift
   \param[in] df_2,                     Speed-dependent line pressure shift
   \param[in] sigma,                    Doppler broadening
   \param[in] f_VC,                     Collisional frequency limit
   \param[in] LM_DF,                    Line mixing frequency shift
   \param[in] LM_Y,                     Line mixing dispersion dependency
   \param[in] LM_G,                     Line mixing added attenuation
   \param[in] cut,                      Is cutoff applied?
   \param[in] calc_phase,               Is dispersion calculated?
   \param[in] calc_src,                 Is source calculated?
 
   \author Stefan Buehler and Axel von Engeln
   \date   2001-01-11 
   
   Changed from pseudo cross sections to true cross sections
   
   \author Stefan Buehler 
   \date   2007-08-08 
   
   Adapted to new Perrin line parameters, treating broadening by different 
   gases explicitly
   
   \author Stefan Buehler
   \date   2012-09-03
   
   Removed LineRecord dependency.
   
   \author Richard Larsson
   \date   2014-10-29
   
*/
void xsec_single_line(  // Output:
    VectorView xsec_accum_attenuation,
    VectorView xsec_accum_source,
    VectorView xsec_accum_phase,
    // Helper variables
    Vector& attenuation,
    Vector& phase,
    Vector& dFa_dx,
    Vector& dFb_dx,
    Vector& dFa_dy,
    Vector& dFb_dy,
    Range& this_f_range,
    Vector& fac,
    Vector& aux,
    // Frequency grid:
    Vector& f_local,
    const ConstVectorView f_grid,
    const Index nf,
    const Numeric cutoff,
    Numeric F0,
    // Line strength:
    Numeric intensity,
    const Numeric part_fct_ratio,
    const Numeric boltzmann_ratio,
    const Numeric abs_nlte_ratio,
    const Numeric src_nlte_ratio,
    const Numeric Isotopologue_Ratio,
    // Atmospheric state
    const Numeric temperature,
    // Line shape:
    const Index ind_ls,
    const Index ind_lsn,
    // Line broadening:
    const Numeric gamma_0,
    const Numeric gamma_2,
    const Numeric eta,
    const Numeric df_0,
    const Numeric df_2,
    const Numeric sigma,
    const Numeric f_VC,
    // Line mixing
    const Numeric LM_DF,
    const Numeric LM_Y,
    const Numeric LM_G,
    // Feature flags
    const bool calc_cut,
    const bool calc_phase,
    const bool calc_partials,
    const bool calc_src) {
  // Copy f_grid to the beginning of f_local. There is one
  // element left at the end of f_local.
  // THIS HAS TO BE INSIDE THE LINE LOOP, BECAUSE THE CUTOFF
  // FREQUENCY IS ALWAYS PUT IN A DIFFERENT PLACE!

  // This will hold the actual number of frequencies to add to
  // xsec later on:
  Index nfl = nf;

  // This will hold the actual number of frequencies for the
  // call to the lineshape functions later on:
  Index nfls = nf;

  // intensity at temperature
  // (calculate the line intensity according to the standard
  // expression given in HITRAN)
  intensity *= part_fct_ratio * boltzmann_ratio;

  // Apply pressure shift:
  //
  // 2017-10-15: Stuart Fox had reported the bug that for mirror
  // lines (negative F0) the shift was in the wrong direction. The
  // if statement below makes sure that the shift goes in the other
  // direction for negative frequency mirror lines, so that they
  // really end up at the negative of the frequency of the original
  // line.
  if (F0 >= 0)
    F0 += (df_0 + LM_DF);
  else
    F0 -= (df_0 + LM_DF);
  // FIXME:  This is probably not compatible with HTP line shape

  // Indices pointing at begin/end frequencies of f_grid or at
  // the elements that have to be calculated in case of cutoff
  Index i_f_min = 0;
  Index i_f_max = nf - 1;

  // cutoff ?
  if (calc_cut) {
    // Check whether we have elements in ls that can be
    // ignored at lower frequencies of f_grid.
    //
    // Loop through all frequencies, finding min value and
    // set all values to zero on that way.
    while (i_f_min < nf && (F0 - cutoff) > f_grid[i_f_min]) {
      //              ls[i_f_min] = 0;
      ++i_f_min;
    }

    // Check whether we have elements in ls that can be
    // ignored at higher frequencies of f_grid.
    //
    // Loop through all frequencies, finding max value and
    // set all values to zero on that way.
    while (i_f_max >= 0 && (F0 + cutoff) < f_grid[i_f_max]) {
      //              ls[i_f_max] = 0;
      --i_f_max;
    }

    // Append the cutoff frequency to f_local:
    ++i_f_max;
    f_local[i_f_max] = F0 + cutoff;

    // Number of frequencies to calculate:
    nfls = i_f_max - i_f_min + 1;  // Add one because indices
                                   // are pointing to first and
                                   // last valid element. This
                                   // is for the lineshape
                                   // calls.
    nfl = nfls - 1;                // This is for xsec.
  } else {
    // Nothing to do here. Note that nfl and nfls are both still set to nf.
  }

  //          cout << "nf, nfl, nfls = " << nf << ", " << nfl << ", " << nfls << ".\n";

  // Maybe there are no frequencies left to compute?  Note that
  // the number that counts here is nfl, since only these are
  // the `real' frequencies, for which xsec is changed. nfls
  // will always be at least one, because it contains the cutoff.
  if (nfl > 0) {
    this_f_range = Range(i_f_min, nfl);
    Range that_f_range(i_f_min, nfls);

    // Calculate the line shape:
    global_data::lineshape_data[ind_ls].Function()(
        attenuation,
        phase,
        dFa_dx,
        dFb_dx,
        dFa_dy,
        dFb_dy,  //partial derivatives
        aux,
        F0,
        gamma_0,
        gamma_2,
        eta,
        0.0,  // FIXME: H-T ls...
        df_2,
        sigma,
        f_VC,
        f_local[that_f_range],
        calc_phase,
        calc_partials);

    // Calculate the chosen normalization factor:
    global_data::lineshape_norm_data[ind_lsn].Function()(
        fac, F0, f_local[that_f_range], temperature);

    // Reset f_local for next loop through now that cutoff is calculated
    if (i_f_max < nf) f_local[i_f_max] = f_grid[i_f_max];

    // Attenuation output
    VectorView this_xsec_attenuation = xsec_accum_attenuation[this_f_range];

    // Optional outputs
    Vector empty_vector;
    VectorView this_xsec_source =
        calc_src ? xsec_accum_source[this_f_range] : empty_vector;
    VectorView this_xsec_phase =
        calc_phase ? xsec_accum_phase[this_f_range] : empty_vector;

    // Store cutoff values
    const Numeric cutoff_attenuation = calc_cut ? attenuation[nfls - 1] : 0.0,
                  cutoff_phase =
                      calc_phase ? calc_cut ? phase[nfls - 1] : 0.0 : 0.0,
                  cutoff_dFa_dx =
                      calc_partials ? calc_cut ? dFa_dx[nfls - 1] : 0.0 : 0.0,
                  cutoff_dFb_dx = calc_phase
                                      ? calc_partials
                                            ? calc_cut ? dFb_dx[nfls - 1] : 0.0
                                            : 0.0
                                      : 0.0,
                  cutoff_dFa_dy =
                      calc_partials ? calc_cut ? dFa_dy[nfls - 1] : 0.0 : 0.0,
                  cutoff_dFb_dy = calc_phase
                                      ? calc_partials
                                            ? calc_cut ? dFb_dy[nfls - 1] : 0.0
                                            : 0.0
                                      : 0.0;

    // If one of these is non-zero, then there are line mixing calculations required
    const bool calc_LM = LM_G != 0 || LM_Y != 0;

    // Line strength is going to be stored in attenuation to ensure the right numeric value.
    const Numeric str_line =
        Isotopologue_Ratio * intensity * (calc_src ? abs_nlte_ratio : 1.0);

    // To make this readable, we loop over the frequencies
    for (Index jj = 0; jj < this_xsec_attenuation.nelem(); jj++) {
      const Numeric str_scale = fac[jj] * str_line;
      if (jj <= nfl) {
        // This ensures that output attenuation is the physically correct attenuation
        if (calc_cut) attenuation[jj] -= cutoff_attenuation;
        attenuation[jj] *= str_scale;

        // This ensures that output phase is the physically correct phase
        if (calc_phase || calc_LM) {
          if (calc_cut) phase[jj] -= cutoff_phase;
          phase[jj] *= str_scale;
        }

        if (calc_partials) {
          if (calc_cut) {
            dFa_dx[jj] -= cutoff_dFa_dx;
            dFa_dy[jj] -= cutoff_dFa_dy;
            if (calc_phase) {
              dFb_dx[jj] -= cutoff_dFb_dx;
              dFb_dy[jj] -= cutoff_dFb_dy;
            }
          }
          dFa_dx[jj] *= str_scale;
          dFa_dy[jj] *= str_scale;
          if (calc_phase) {
            dFb_dx[jj] *= str_scale;
            dFb_dy[jj] *= str_scale;
          }
        }

        // Attenuation is by tradition added to total attenuation from here
        if (calc_LM && calc_src) {
          const Numeric tmp = (1. + LM_G) * attenuation[jj] + LM_Y * phase[jj];
          this_xsec_source[jj] += (src_nlte_ratio / abs_nlte_ratio - 1.0) * tmp;
          this_xsec_attenuation[jj] += tmp;
        } else if (calc_src) {
          this_xsec_source[jj] +=
              (src_nlte_ratio / abs_nlte_ratio - 1.0) * attenuation[jj];
          this_xsec_attenuation[jj] += attenuation[jj];
        } else if (calc_LM)
          this_xsec_attenuation[jj] +=
              (1. + LM_G) * attenuation[jj] + LM_Y * phase[jj];
        else
          this_xsec_attenuation[jj] += attenuation[jj];

        // Phase follows the attenuation practice
        if (calc_phase && calc_LM)
          this_xsec_phase[jj] +=
              (1. + LM_G) * phase[jj] - LM_Y * attenuation[jj];
        else if (calc_phase)
          this_xsec_phase[jj] += phase[jj];

        // Partials (Assuming F(v,v0) = C*F'(v,v0), where F' is from our lineshape functions, then the below is
        // only the C*dF'/dt part.  The dC/dt*F' part remains to be calculated.  The CF is output as
        // attenuation and phase variables [NOTE: attenuation and phase are without line mixing corrections!])
        if (calc_partials && calc_LM) {
          const Numeric orig_dFa_dx = dFa_dx[jj];
          const Numeric orig_dFa_dy = dFa_dy[jj];

          dFa_dx[jj] = (1. + LM_G) * orig_dFa_dx + LM_Y * dFb_dx[jj];
          dFb_dx[jj] = (1. + LM_G) * dFb_dx[jj] - LM_Y * orig_dFa_dx;
          dFa_dy[jj] = (1. + LM_G) * orig_dFa_dy + LM_Y * dFb_dy[jj];
          dFb_dy[jj] = (1. + LM_G) * dFb_dy[jj] - LM_Y * orig_dFa_dy;
        }
      }
    }
  }
}

/** A little helper function to convert energy from units of
    wavenumber (cm^-1) to Joule (J). 

    This is used when reading HITRAN or JPL catalogue files, which
    have the lower state energy in cm^-1.

    \return Energy in J.
    \param[in]  e Energy in cm^-1.

    \author Stefan Buehler
    \date   2001-06-26 */
Numeric wavenumber_to_joule(Numeric e) {
  // Planck constant [Js]
  extern const Numeric PLANCK_CONST;

  // Speed of light [m/s]
  extern const Numeric SPEED_OF_LIGHT;

  // Constant to convert lower state energy from cm^-1 to J
  const Numeric lower_energy_const = PLANCK_CONST * SPEED_OF_LIGHT * 1E2;

  return e * lower_energy_const;
}

//======================================================================
//             Functions related to species
//======================================================================

//! Return species index for given species name.
/*! 
  This is useful in connection with other functions that need a species
  index.

  \see find_first_species_tg.

  \param[in] name Species name.

  \return Species index, -1 means not found.

  \author Stefan Buehler
  \date   2003-01-13
*/
Index species_index_from_species_name(String name) {
  using global_data::SpeciesMap;

  // For the return value:
  Index mspecies;

  // Trim whitespace
  name.trim();

  //  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<String, Index>::const_iterator mi = SpeciesMap.find(name);
  if (mi != SpeciesMap.end()) {
    // Ok, we've found the species. Set mspecies.
    mspecies = mi->second;
  } else {
    // The species does not exist!
    mspecies = -1;
  }

  return mspecies;
}

//! Return species name for given species index.
/*!
 This is useful in connection with other functions that use a species
 index.
 
 Does an assertion that the index really corresponds to a species.
 
 \param[in] spec_ind Species index.
 
 \return Species name
 
 \author Stefan Buehler
 \date   2013-01-04
 */
String species_name_from_species_index(const Index spec_ind) {
  // Species lookup data:
  using global_data::species_data;

  // Assert that spec_ind is inside species data. (This is an assertion,
  // because species indices should never be user input, but set by the
  // program automatically, based on species names.)
  assert(spec_ind >= 0);
  assert(spec_ind < species_data.nelem());

  // A reference to the relevant record of the species data:
  const SpeciesRecord& spr = species_data[spec_ind];

  return spr.Name();
}

//======================================================================
//        Functions to convert the accuracy index to ARTS units
//======================================================================

// ********* for HITRAN database *************
//convert HITRAN index for line position accuracy to ARTS
//units (Hz).

void convHitranIERF(Numeric& mdf, const Index& df) {
  switch (df) {
    case 0: {
      mdf = -1;
      break;
    }
    case 1: {
      mdf = 30 * 1E9;
      break;
    }
    case 2: {
      mdf = 3 * 1E9;
      break;
    }
    case 3: {
      mdf = 300 * 1E6;
      break;
    }
    case 4: {
      mdf = 30 * 1E6;
      break;
    }
    case 5: {
      mdf = 3 * 1E6;
      break;
    }
    case 6: {
      mdf = 0.3 * 1E6;
      break;
    }
  }
}

//convert HITRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convHitranIERSH(Numeric& mdh, const Index& dh) {
  switch (dh) {
    case 0: {
      mdh = -1;
      break;
    }
    case 1: {
      mdh = -1;
      break;
    }
    case 2: {
      mdh = -1;
      break;
    }
    case 3: {
      mdh = 30;
      break;
    }
    case 4: {
      mdh = 20;
      break;
    }
    case 5: {
      mdh = 10;
      break;
    }
    case 6: {
      mdh = 5;
      break;
    }
    case 7: {
      mdh = 2;
      break;
    }
    case 8: {
      mdh = 1;
      break;
    }
  }
  mdh = mdh / 100;
}

// ********* for MYTRAN database *************
//convert MYTRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convMytranIER(Numeric& mdh, const Index& dh) {
  switch (dh) {
    case 0: {
      mdh = 200;
      break;
    }
    case 1: {
      mdh = 100;
      break;
    }
    case 2: {
      mdh = 50;
      break;
    }
    case 3: {
      mdh = 30;
      break;
    }
    case 4: {
      mdh = 20;
      break;
    }
    case 5: {
      mdh = 10;
      break;
    }
    case 6: {
      mdh = 5;
      break;
    }
    case 7: {
      mdh = 2;
      break;
    }
    case 8: {
      mdh = 1;
      break;
    }
    case 9: {
      mdh = 0.5;
      break;
    }
  }
  mdh = mdh / 100;
}

ostream& operator<<(ostream& os, const LineshapeSpec& lsspec) {
  os << "LineshapeSpec Index: " << lsspec.Ind_ls()
     << ", NormIndex: " << lsspec.Ind_lsn() << ", Cutoff: " << lsspec.Cutoff()
     << endl;

  return os;
}

/**
 *  
 *  This will work as a wrapper for linemixing when abs_species contain relevant data.
 *  The funciton will only pass on arguments to xsec_species if there is no linemixing.
 *  
 *  \param[out] xsec_attenuation    Cross section of one tag group. This is now the
 *                              true attenuation cross section in units of m^2.
 *  \param[out] xsec_source         Cross section of one tag group. This is now the
 *                              true source cross section in units of m^2.
 *  \param[out] xsec_phase          Cross section of one tag group. This is now the
 *                              true phase cross section in units of m^2.
 *  \param[in] f_grid               Frequency grid.
 *  \param[in] abs_p                Pressure grid.
 *  \param[in] abs_t                Temperatures associated with abs_p.
 *  \param[in] all_vmrs             Gas volume mixing ratios [nspecies, np].
 *  \param[in] abs_species          Species tags for all species.
 *  \param[in] this_species         Index of the current species in abs_species.
 *  \param[in] abs_lines            The spectroscopic line list.
 *  \param[in] ind_ls               Index to lineshape function.
 *  \param[in] ind_lsn              Index to lineshape norm.
 *  \param[in] cutoff               Lineshape cutoff.
 *  \param[in] isotopologue_ratios  Isotopologue ratios.
 * 
 *  \author Richard Larsson
 *  \date   2013-04-24
 * 
 */
void xsec_species_line_mixing_wrapper(
    MatrixView xsec_attenuation,
    MatrixView xsec_source,
    MatrixView xsec_phase,
    ArrayOfMatrix& partial_xsec_attenuation,
    ArrayOfMatrix& partial_xsec_source,
    ArrayOfMatrix& partial_xsec_phase,
    const ArrayOfRetrievalQuantity& flag_partials,
    const ArrayOfIndex& flag_partials_position,
    ConstVectorView f_grid,
    ConstVectorView abs_p,
    ConstVectorView abs_t,
    ConstMatrixView abs_t_nlte,
    ConstMatrixView all_vmrs,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Index this_species,
    const ArrayOfLineRecord& abs_lines,
    const Numeric H_magntitude_Zeeman,
    const Index ind_ls,
    const Index ind_lsn,
    const Numeric cutoff,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions) {
  // Optional paths through the code...
  const bool cut = (cutoff != -1) ? true : false;
  const bool calc_partials = flag_partials_position.nelem();
  const bool calc_partials_phase =
      partial_xsec_phase.nelem() == partial_xsec_attenuation.nelem();
  const bool do_zeeman = false;
  const bool do_lm =
      abs_species[this_species][0].LineMixing() != SpeciesTag::LINE_MIXING_OFF;
  const bool no_ls_partials = supports_LBL_without_phase(flag_partials);

  const bool we_need_phase = do_zeeman || do_lm || calc_partials_phase;
  const bool we_need_partials = calc_partials && !no_ls_partials;

  // Must have the phase for two of the options
  using global_data::lineshape_data;
  if (!lineshape_data[ind_ls].Phase() && we_need_phase) {
    std::ostringstream os;
    os << "This is an error message. You are using "
       << lineshape_data[ind_ls].Name() << ".\n"
       << "This line shape does not include phase in its calculations and\nis therefore invalid for "
       << "line mixing, Zeeman effect, and certain partial derivatives.\nYou are using one of these or you should not see this error.\n";
    throw std::runtime_error(os.str());
  }

  if (we_need_partials && !lineshape_data[ind_ls].Partials()) {
    std::ostringstream os;
    os << "This is an error message. You are using "
       << lineshape_data[ind_ls].Name() << ".\n"
       << "Your selected *jacobian_quantities* requires that the line shape returns partial\n"
       << "derivatives.";

    throw std::runtime_error(os.str());
  }

  extern const Numeric DOPPLER_CONST;
  static const Numeric doppler_const = DOPPLER_CONST;

  // Check that the frequency grid is sorted in the case of lineshape
  // with cutoff. Duplicate frequency values are allowed.
  if (cut) {
    if (!is_sorted(f_grid)) {
      std::ostringstream os;
      os << "If you use a lineshape function with cutoff, your\n"
         << "frequency grid *f_grid* must be sorted.\n"
         << "(Duplicate values are allowed.)";
      throw std::runtime_error(os.str());
    }
  }

  // Check that all temperatures are non-negative
  bool negative = false;

  for (Index i = 0; !negative && i < abs_t.nelem(); i++) {
    if (abs_t[i] < 0.) negative = true;
  }

  if (negative) {
    std::ostringstream os;
    os << "abs_t contains at least one negative temperature value.\n"
       << "This is not allowed.";
    throw std::runtime_error(os.str());
  }

  if (abs_t.nelem() != abs_p.nelem()) {
    std::ostringstream os;
    os << "Variable abs_t must have the same dimension as abs_p.\n"
       << "abs_t.nelem() = " << abs_t.nelem() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw std::runtime_error(os.str());
  }

  // all_vmrs should have dimensions [nspecies, np]:

  if (all_vmrs.ncols() != abs_p.nelem()) {
    std::ostringstream os;
    os << "Number of columns of all_vmrs must match abs_p.\n"
       << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw std::runtime_error(os.str());
  }

  // Check that the dimension of xsec is indeed [f_grid.nelem(),
  // abs_p.nelem()]:
  if (xsec_attenuation.nrows() != f_grid.nelem() ||
      xsec_attenuation.ncols() != abs_p.nelem()) {
    std::ostringstream os;
    os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
       << "[xsec_attenuation.nrows(),xsec_attenuation.ncols()] = ["
       << xsec_attenuation.nrows() << ", " << xsec_attenuation.ncols() << "]\n"
       << "f_grid.nelem() = " << f_grid.nelem() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw std::runtime_error(os.str());
  }

  // Check that the dimension of xsec source is [f_grid.nelem(), abs_p.nelem()]:
  bool calc_src;
  if (xsec_source.nrows() != f_grid.nelem() ||
      xsec_source.ncols() != abs_p.nelem()) {
    if (xsec_source.nrows() != 0 || xsec_source.ncols() != 0) {
      std::ostringstream os;
      os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()] or [0,0].\n"
         << "[xsec_source.nrows(),xsec_source.ncols()] = ["
         << xsec_source.nrows() << ", " << xsec_source.ncols() << "]\n"
         << "f_grid.nelem() = " << f_grid.nelem() << '\n'
         << "abs_p.nelem() = " << abs_p.nelem();
      throw std::runtime_error(os.str());
    } else {
      calc_src = 0;
    }
  } else {
    calc_src = 1;
  }

  if (xsec_phase.nrows() != f_grid.nelem() ||
      xsec_phase.ncols() != abs_p.nelem()) {
    std::ostringstream os;
    os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
       << "[xsec_phase.nrows(),xsec_phase.ncols()] = [" << xsec_phase.nrows()
       << ", " << xsec_phase.ncols() << "]\n"
       << "f_grid.nelem() = " << f_grid.nelem() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw std::runtime_error(os.str());
  }

  //Helper var
  Vector attenuation(f_grid.nelem() + 1), phase(f_grid.nelem() + 1),
      aux(f_grid.nelem() + 1), fac(f_grid.nelem() + 1),
      f_local(f_grid.nelem() + 1);
  f_local[Range(0, f_grid.nelem())] = f_grid;

  // Set this to zero in the case of VoigtKuntz6 or other lines hapes that work in special circumstances
  if ((not we_need_phase) and we_need_partials) phase = 0.0;

#pragma omp parallel for if (!arts_omp_in_parallel()) \
    firstprivate(attenuation, phase, fac, f_local, aux)
  for (Index jj = 0; jj < abs_p.nelem(); jj++) {
    Vector empty_vector;
    const Numeric& t = abs_t[jj];
    const Numeric& p = abs_p[jj];
    ConstVectorView t_nlte = calc_src ? abs_t_nlte(joker, jj) : empty_vector;
    ConstVectorView vmrs = all_vmrs(joker, jj);

    // Setup for calculating the partial derivatives
    Vector dFa_dx, dFb_dx, dFa_dy, dFb_dy;
    if (calc_partials &&
        !no_ls_partials)  //Only size them if partials are wanted
    {
      dFa_dx.resize(f_grid.nelem() + 1);
      dFb_dx.resize(f_grid.nelem() + 1);
      dFa_dy.resize(f_grid.nelem() + 1);
      dFb_dy.resize(f_grid.nelem() + 1);
    }

    // Simple caching of partition function to avoid recalculating things.
    Numeric qt_cache = -1, qref_cache = -1;
    Index iso_cache = -1;
    Numeric line_t_cache = -1, dqt_dT_cache = -1.;

    for (Index ii = 0; ii < abs_lines.nelem(); ii++) {
      // These lines should be ignored by user request
      if (abs_lines[ii].LineMixingByBand()) continue;

      if (calc_src) {
        if (abs_lines[ii].GetLinePopulationType() not_eq
                LinePopulationType::ByLTE and
            abs_lines[ii].GetLinePopulationType() not_eq
                LinePopulationType::ByVibrationalTemperatures) {
          throw std::runtime_error(
              "Bad data seen in xsec_species_line_mixing_wrapper.  Please use more modern function.");
        }
      }

      // Pressure broadening parameters
      // Prepare pressure broadening parameters
      auto X = abs_lines[ii].GetShapeParams(t, p, vmrs, abs_species);
      const Numeric& G0 = X.G0;
      const Numeric& D0 = X.D0;
      const Numeric& G2 = X.G2;
      const Numeric& D2 = X.D2;
      const Numeric& FVC = X.FVC;
      const Numeric& ETA = X.ETA;
      const Numeric& Y = X.Y;
      const Numeric& G = X.G;
      const Numeric& DV = X.DV;

      // Check the cache is the temperature of the line and the isotope is the same to avoid recalculating the partition sum
      if (iso_cache != abs_lines[ii].Isotopologue() ||
          line_t_cache != abs_lines[ii].Ti0()) {
        iso_cache = abs_lines[ii].Isotopologue();
        line_t_cache = abs_lines[ii].Ti0();

        qt_cache = -1;
        qref_cache = -1;
        dqt_dT_cache = -1.;
      }

      // Prepare line strength scaling
      Numeric partition_ratio, K1, K2, abs_nlte_ratio, src_nlte_ratio;
      GetLineScalingData(
          qt_cache,
          qref_cache,
          partition_ratio,
          K1,
          K2,
          abs_nlte_ratio,
          src_nlte_ratio,
          partition_functions.getParamType(abs_lines[ii].Species(),
                                           abs_lines[ii].Isotopologue()),
          partition_functions.getParam(abs_lines[ii].Species(),
                                       abs_lines[ii].Isotopologue()),
          t,
          abs_lines[ii].Ti0(),
          abs_lines[ii].F(),
          abs_lines[ii].Elow(),
          calc_src,
          abs_lines[ii].Evlow(),
          abs_lines[ii].Evupp(),
          abs_lines[ii].NLTELowerIndex(),
          abs_lines[ii].NLTEUpperIndex(),
          t_nlte);

      // Doppler broadening
      const Numeric sigma = abs_lines[ii].F() * doppler_const *
                            sqrt(t / abs_lines[ii].IsotopologueData().Mass());

      Range this_f_range(0, 0);

      xsec_single_line(  // OUTPUT
          xsec_attenuation(joker, jj),
          calc_src ? xsec_source(joker, jj) : empty_vector,
          xsec_phase(joker, jj),
          // HELPER
          attenuation,
          phase,
          dFa_dx,
          dFb_dx,
          dFa_dy,
          dFb_dy,
          this_f_range,
          fac,
          aux,
          // FREQUENCY
          f_local,
          f_grid,
          f_grid.nelem(),
          cutoff,
          abs_lines[ii].F(),
          // LINE STRENGTH
          abs_lines[ii].I0(),
          partition_ratio,
          K1 * K2,
          abs_nlte_ratio,
          src_nlte_ratio,
          isotopologue_ratios
              .getParam(abs_lines[ii].Species(),
                        abs_lines[ii].Isotopologue())[0]
              .data[0],
          // ATMOSPHERIC TEMPERATURE
          t,
          // LINE SHAPE
          ind_ls,
          ind_lsn,
          // LINE BROADENING
          G0,
          G2,
          ETA,
          D0,
          D2,
          sigma,
          FVC,
          // LINE MIXING
          DV,
          Y,
          G,
          // FEATURE FLAGS
          cutoff > 0,
          we_need_phase,
          calc_partials && !no_ls_partials,
          calc_src);

      if (calc_partials) {
        // These needs to be calculated and returned when Temperature is in list
        Numeric dG0dT, dD0dT, dYdT, dGdT, dDVdT, dQ_dT, dK2_dT,
            dabs_nlte_ratio_dT = 0.0, atm_tv_low, atm_tv_upp;
        if (do_temperature_jacobian(flag_partials)) {
          auto dX = abs_lines[ii].GetShapeParams_dT(t, p, vmrs, abs_species);
          dG0dT = dX.G0;
          dD0dT = dX.D0;
          dYdT = dX.Y;
          dGdT = dX.G;
          dDVdT = dX.DV;

          GetLineScalingData_dT(
              dqt_dT_cache,
              dK2_dT,
              dQ_dT,
              dabs_nlte_ratio_dT,
              atm_tv_low,
              atm_tv_upp,
              qt_cache,
              abs_nlte_ratio,
              partition_functions.getParamType(abs_lines[ii].Species(),
                                               abs_lines[ii].Isotopologue()),
              partition_functions.getParam(abs_lines[ii].Species(),
                                           abs_lines[ii].Isotopologue()),
              t,
              abs_lines[ii].Ti0(),
              temperature_perturbation(flag_partials),
              abs_lines[ii].F(),
              calc_src,
              abs_lines[ii].Evlow(),
              abs_lines[ii].Evupp(),
              abs_lines[ii].NLTELowerIndex(),
              abs_lines[ii].NLTEUpperIndex(),
              t_nlte);
        }
        
        for (Index iq = 0; iq < flag_partials_position.nelem(); iq++) {
          if (is_line_parameter(flag_partials[flag_partials_position[iq]]))
            throw std::runtime_error(
                "Sorry, there has been a hard deprecation of old xsec line param derivatives\nPlease use a modern xsec method");
        }

        // Gather all new partial derivative calculations in this function
        partial_derivatives_lineshape_dependency(
            partial_xsec_attenuation,
            partial_xsec_phase,
            partial_xsec_source,
            flag_partials,
            flag_partials_position,
            attenuation,
            phase,
            fac,
            dFa_dx,
            dFb_dx,
            dFa_dy,
            dFb_dy,
            f_grid,
            this_f_range,
            //Temperature
            t,
            sigma,
            K2,
            dK2_dT,
            abs_nlte_ratio,  //K3
            dabs_nlte_ratio_dT,
            src_nlte_ratio,  //K4
            // Line parameters
            abs_lines[ii].F(),
            abs_lines[ii].I0(),
            abs_lines[ii].Ti0(),
            abs_lines[ii].Elow(),
            abs_lines[ii].Evlow(),
            abs_lines[ii].Evupp(),
            (abs_lines[ii].NLTELowerIndex() > -1 and calc_src)
                ? t_nlte[abs_lines[ii].NLTELowerIndex()]
                : -1.0,
            (abs_lines[ii].NLTEUpperIndex() > -1 and calc_src)
                ? t_nlte[abs_lines[ii].NLTEUpperIndex()]
                : -1.0,
            Y,
            dYdT,
            G,
            dGdT,
            DV,
            dDVdT,
            abs_lines[ii].QuantumIdentity(),
            // LINE SHAPE
            ind_ls,
            ind_lsn,
            D0,  //FIXME: For H-T, this part is difficult
            dD0dT,
            G0,
            dG0dT,
            // Partition data parameters
            dQ_dT,
            // Magnetic variables
            0.0,
            H_magntitude_Zeeman,
            false,
            // Programming
            jj,
            calc_partials_phase,
            calc_src);
      }
    }
  }
}

/*! cross-section replacement computer 
 *  
 * This will work as the interface for all line-by-line computations 
 * lacking special demands
 *  
 *  \param[in,out] xsec         Cross section of one tag group. This is now the
 *                              true attenuation cross section in units of m^2.
 *  \param[in,out] source       Cross section of one tag group. This is now the
 *                              true source cross section in units of m^2.
 *  \param[in,out] phase        Cross section of one tag group. This is now the
 *                              true phase cross section in units of m^2.
 *  \param[in,out] dxsec        Partial derivatives of xsec.
 *  \param[in,out] dsource      Partial derivatives of source.
 *  \param[in,out] dphase       Partial derivatives of phase.
 * 
 *  \param[in] flag_partials        Partial derivatives flags.
 *  \param[in] f_grid               Frequency grid.
 *  \param[in] abs_p                Pressure grid.
 *  \param[in] abs_t                Temperatures associated with abs_p.
 *  \param[in] abs_t_nlte           Non-lte temperatures for various energy levels.
 *  \param[in] all_vmrs             Gas volume mixing ratios [nspecies, np].
 *  \param[in] abs_species          Species tags for all species.
 *  \param[in] this_species         Index of the current species in abs_species.
 *  \param[in] abs_lines            The spectroscopic line list.
 *  \param[in] Z_DF                 The Zeeman line center shift over the magnitude of the magnetic field.
 *  \param[in] H_magntitude_Zeeman  The magnitude of the magnetic field required by Zeeman effect.
 *  \param[in] lm_p_lim             Line mixing pressure limit
 *  \param[in] isotopologue_ratios  Isotopologue ratios.
 *  \param[in] partition_functions  Partition functions.
 *  \param[in] verbosity            Verbosity level.
 * 
 *  \author Richard Larsson
 *  \date   2013-04-24
 * 
 */
void xsec_species2(Matrix& xsec,
                   Matrix& source,
                   Matrix& phase,
                   ArrayOfMatrix& dxsec_dx,
                   ArrayOfMatrix& dsource_dx,
                   ArrayOfMatrix& dphase_dx,
                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                   const ArrayOfIndex& jacobian_propmat_positions,
                   const Vector& f_grid,
                   const Vector& abs_p,
                   const Vector& abs_t,
                   const Matrix& abs_nlte,
                   const Matrix& all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const ArrayOfLineRecord& abs_lines,
                   const SpeciesAuxData& isotopologue_ratios,
                   const SpeciesAuxData& partition_functions) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = abs_lines.nelem();  // number of lines in the catalog
  const Index nj =
      jacobian_propmat_positions.nelem();  // number of partial derivatives
  const Index nt = source.nrows();         // number of energy levels in NLTE

  // Type of problem
  const bool do_nonlte = nt;
  const bool do_jacobi = nj;
  const bool do_temperature = do_temperature_jacobian(jacobian_quantities);

  // Test if the size of the problem is 0
  if (not np or not nf or not nl) return;

  // Move the problem to Eigen-library types
  const auto f_grid_eigen = MapToEigen(f_grid);

  // Parallelization data holder (skip if in parallel or there are too few threads)
  const Index nthread =
      arts_omp_in_parallel()
          ? 1
          : arts_omp_get_max_threads() < nl ? arts_omp_get_max_threads() : 1;
  std::vector<Eigen::VectorXcd> Fsum(nthread, Eigen::VectorXcd(nf));
  std::vector<Eigen::MatrixXcd> dFsum(nthread, Eigen::MatrixXcd(nf, nj));
  std::vector<Eigen::VectorXcd> Nsum(nthread, Eigen::VectorXcd(nf));
  std::vector<Eigen::MatrixXcd> dNsum(nthread, Eigen::MatrixXcd(nf, nj));

  for (Index ip = 0; ip < np; ip++) {
    // Constants for this level
    const Numeric& temperature = abs_t[ip];
    const Numeric& pressure = abs_p[ip];

    // Quasi-constants for this level, defined here to speed up later computations
    Index this_iso = -1;  // line isotopologue number
    Numeric t0 = -1;      // line temperature
    Numeric dc = 0, ddc_dT = 0, qt = 0, qt0 = 1,
            dqt_dT = 0;  // Doppler and partition functions

    // Line shape constants
    LineShape::Model line_shape_default_model;
    LineShape::Model& line_shape_model{line_shape_default_model};
    Vector line_shape_vmr(0);

    // Reset sum-operators
    for (Index ithread = 0; ithread < nthread; ithread++) {
      Fsum[ithread].setZero();
      dFsum[ithread].setZero();
      Nsum[ithread].setZero();
      dNsum[ithread].setZero();
    }

#pragma omp parallel for if (nthread > 1) schedule(guided, 4) \
    firstprivate(this_iso,                                    \
                 t0,                                          \
                 dc,                                          \
                 ddc_dT,                                      \
                 qt,                                          \
                 qt0,                                         \
                 dqt_dT,                                      \
                 line_shape_model,                            \
                 line_shape_vmr)
    for (Index il = 0; il < abs_lines.nelem(); il++) {
      const auto& line = abs_lines[il];

      // Local compute variables
      thread_local Eigen::VectorXcd F(nf);
      thread_local Eigen::MatrixXcd dF(nf, nj);
      thread_local Eigen::VectorXcd N(nf);
      thread_local Eigen::MatrixXcd dN(nf, nj);
      thread_local Eigen::
          Matrix<Complex, Eigen::Dynamic, Linefunctions::ExpectedDataSize()>
              data(nf, Linefunctions::ExpectedDataSize());
      thread_local Index start, nelem;

      // Partition function depends on isotopologue and line temperatures.
      // Both are commonly constant in a single catalog.  They are, however,
      // allowed to change so we must check that they do not
      if (line.Isotopologue() not_eq this_iso or line.Ti0() not_eq t0) {
        t0 = line.Ti0();

        partition_function(
            qt0,
            qt,
            t0,
            temperature,
            partition_functions.getParamType(line.Species(),
                                             line.Isotopologue()),
            partition_functions.getParam(line.Species(), line.Isotopologue()));

        if (do_temperature)
          dpartition_function_dT(dqt_dT,
                                 qt,
                                 temperature,
                                 temperature_perturbation(jacobian_quantities),
                                 partition_functions.getParamType(
                                     line.Species(), line.Isotopologue()),
                                 partition_functions.getParam(
                                     line.Species(), line.Isotopologue()));

        if (line.Isotopologue() not_eq this_iso) {
          this_iso = line.Isotopologue();
          dc = Linefunctions::DopplerConstant(temperature,
                                              line.IsotopologueData().Mass());
          if (do_temperature)
            ddc_dT = Linefunctions::dDopplerConstant_dT(temperature, dc);
        }
      }

      if (not line_shape_model.same_broadening_species(
              line.GetLineShapeModel())) {
        line_shape_model = line.GetLineShapeModel();
        line_shape_vmr = line_shape_model.vmrs(
            all_vmrs(joker, ip), abs_species, line.QuantumIdentity());
      }

      Linefunctions::set_cross_section_for_single_line(
          F,
          dF,
          N,
          dN,
          data,
          start,
          nelem,
          f_grid_eigen,
          line,
          jacobian_quantities,
          jacobian_propmat_positions,
          line_shape_vmr,
          nt ? abs_nlte(joker, ip) : Vector(0),
          pressure,
          temperature,
          dc,
          isotopologue_ratios.getParam(line.Species(), this_iso)[0].data[0],
          0.0,
          0.0,
          ddc_dT,
          qt,
          dqt_dT,
          qt0);

      auto ithread = arts_omp_get_thread_num();
      Fsum[ithread].segment(start, nelem).noalias() += F.segment(start, nelem);
      if (do_jacobi)
        dFsum[ithread].middleRows(start, nelem).noalias() +=
            dF.middleRows(start, nelem);
      if (do_nonlte)
        Nsum[ithread].segment(start, nelem).noalias() +=
            N.segment(start, nelem);
      if (do_nonlte and do_jacobi)
        dNsum[ithread].middleRows(start, nelem).noalias() +=
            dN.middleRows(start, nelem);
    }

    // Sum all the threaded results
    for (Index ithread = 0; ithread < nthread; ithread++) {
      // absorption cross-section
      MapToEigen(xsec).col(ip).noalias() += Fsum[ithread].real();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dxsec_dx[j]).col(ip).noalias() +=
            dFsum[ithread].col(j).real();

      // phase cross-section
      if (not phase.empty()) {
        MapToEigen(phase).col(ip).noalias() += Fsum[ithread].imag();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dphase_dx[j]).col(ip).noalias() +=
              dFsum[ithread].col(j).imag();
      }

      // source ratio cross-section
      if (do_nonlte) {
        MapToEigen(source).col(ip).noalias() += Nsum[ithread].real();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dsource_dx[j]).col(ip).noalias() +=
              dNsum[ithread].col(j).real();
      }
    }
  }
}


const SpeciesRecord& SpeciesDataOfLines(const AbsorptionLines& lines)
{
  return global_data::species_data[lines.Species()];
}


/** Cross-section algorithm
 * 
 *  \author Richard Larsson
 *  \date   2019-10-10
 * 
 */
void xsec_species3(Matrix& xsec,
                   Matrix& source,
                   Matrix& phase,
                   ArrayOfMatrix& dxsec_dx,
                   ArrayOfMatrix& dsource_dx,
                   ArrayOfMatrix& dphase_dx,
                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                   const ArrayOfIndex& jacobian_propmat_positions,
                   const Vector& f_grid,
                   const Vector& abs_p,
                   const Vector& abs_t,
                   const EnergyLevelMap& abs_nlte,
                   const Matrix& all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const AbsorptionLines& lines,
                   const Numeric& isot_ratio,
                   const SpeciesAuxData::AuxType& partfun_type,
                   const ArrayOfGriddedField1& partfun_data) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = lines.NumLines();  // number of lines in the catalog
  const Index nj =
      jacobian_propmat_positions.nelem();  // number of partial derivatives
  const Index nt = source.nrows();         // number of energy levels in NLTE

  // Type of problem
  const bool do_nonlte = nt;

  Linefunctions::InternalData scratch(nf, nj), sum(nf, nj);
  
  // Test if the size of the problem is 0
  if (not np or not nf or not nl) return;
  
  // Constant for all lines
  const Numeric QT0 = single_partition_function(lines.T0(), partfun_type, partfun_data);
  
  for (Index ip = 0; ip < np; ip++) {
    // Constants for this level
    const Numeric& temperature = abs_t[ip];
    const Numeric& pressure = abs_p[ip];
    
    // Constants for this level
    const Numeric QT = single_partition_function(temperature, partfun_type, partfun_data);
    const Numeric dQTdT = dsingle_partition_function_dT(QT, temperature, temperature_perturbation(jacobian_quantities), partfun_type, partfun_data);
    const Numeric DC = Linefunctions::DopplerConstant(temperature, lines.SpeciesMass());
    const Numeric dDCdT = Linefunctions::dDopplerConstant_dT(temperature, DC);
    const Vector line_shape_vmr = lines.BroadeningSpeciesVMR(all_vmrs(joker, ip), abs_species);
    
    Linefunctions::set_cross_section_of_band(
      scratch,
      sum,
      f_grid,
      lines,
      jacobian_quantities,
      jacobian_propmat_positions,
      line_shape_vmr,
      abs_nlte[ip],
      pressure,
      temperature,
      isot_ratio,
      0,
      DC,
      dDCdT,
      QT,
      dQTdT,
      QT0,
      true);
    
    // absorption cross-section
    MapToEigen(xsec).col(ip).noalias() += sum.F.real();
    for (Index j = 0; j < nj; j++)
      MapToEigen(dxsec_dx[j]).col(ip).noalias() +=
      sum.dF.col(j).real();
    
    // phase cross-section
    if (not phase.empty()) {
      MapToEigen(phase).col(ip).noalias() += sum.F.imag();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dphase_dx[j]).col(ip).noalias() +=
        sum.dF.col(j).imag();
    }
    
    // source ratio cross-section
    if (do_nonlte) {
      MapToEigen(source).col(ip).noalias() += sum.N.real();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dsource_dx[j]).col(ip).noalias() +=
        sum.dN.col(j).real();
    }
  }
}
