/* Copyright (C) 2019
 Richard Larsson <larsson@mps.mpg.de>

 
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

/** Contains the absorption namespace
 * @file   lineshapemodel.cc
 * @author Richard Larsson
 * @date   2019-09-07
 * 
 * @brief  Contains the absorption lines implementation
 * 
 * This namespace contains classes to deal with absorption lines
 **/

#include "absorptionlines.h"

#include "absorption.h"
#include "global_data.h"
#include "quantum_parser_hitran.h"

Rational Absorption::Lines::LowerQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
  for(size_t i=0; i<mlocalquanta.size(); i++)
    if(mlocalquanta[i] == qnt)
      return mlines[k].LowerQuantumNumber(i);
  
  return mquantumidentity.LowerQuantumNumber(qnt);
}

Rational Absorption::Lines::UpperQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
  for(size_t i=0; i<mlocalquanta.size(); i++)
    if(mlocalquanta[i] == qnt)
      return mlines[k].UpperQuantumNumber(i);
  
  return mquantumidentity.UpperQuantumNumber(qnt);
}

bool Absorption::Lines::InLowerLevel(size_t k,
                                     const QuantumIdentifier& qid) const noexcept {
  if(mquantumidentity.Species() not_eq qid.Species())
    return false;
  else if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
    return false;
  else if(qid.Type() == QuantumIdentifier::ALL)
    return true;
  else if(qid.Type() == QuantumIdentifier::NONE)
    return false;
  else if(qid.Type() == QuantumIdentifier::TRANSITION)
    return false;
  else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL)
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(qid.EnergyLevelQuantumNumber(mlocalquanta[i]).isDefined() and
         qid.EnergyLevelQuantumNumber(mlocalquanta[i]) not_eq mlines[k].LowerQuantumNumber(i))
        return false;
  
  return qid.InLower(mquantumidentity.LowerQuantumId());
}

bool Absorption::Lines::InUpperLevel(size_t k,
                                     const QuantumIdentifier& qid) const noexcept {
  if(mquantumidentity.Species() not_eq qid.Species())
    return false;
  else if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
    return false;
  else if(qid.Type() == QuantumIdentifier::ALL)
    return true;
  else if(qid.Type() == QuantumIdentifier::NONE)
    return false;
  else if(qid.Type() == QuantumIdentifier::TRANSITION)
    return false;
  else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL)
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(qid.EnergyLevelQuantumNumber(mlocalquanta[i]).isDefined() and
         qid.EnergyLevelQuantumNumber(mlocalquanta[i]) not_eq mlines[k].UpperQuantumNumber(i))
        return false;
  
  return qid.InUpper(mquantumidentity.UpperQuantumId());
}

LineShape::Output Absorption::Lines::ShapeParameters(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
  auto x = mlines[k].LineShape().GetParams(T, mT0, P, vmrs);
  
  if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
  
  return x;
}

LineShape::Output Absorption::Lines::ShapeParameters_dT(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
  auto x = mlines[k].LineShape().GetTemperatureDerivs(T, mT0, P, vmrs);
  
  if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
  
  return x;
}

Index Absorption::Lines::LineShapePos(const Index& spec) const noexcept {
  // Is always first if this is self and self broadening exists
  if(mselfbroadening and spec == mquantumidentity.Species())
    return 0;
  
  // First and last might be artificial so they should not be checked
  for(Index i=Index(mselfbroadening); i<Index(mbroadeningspecies.size())-Index(mbathbroadening); i++) {
    if(spec == mbroadeningspecies[i].Species())
      return Index(i);
  }
  
  // At this point, the ID is not explicitly among the broadeners, but bath broadening means its VMR still might matter
  if(mbathbroadening)
    return Index(mbroadeningspecies.size())-Index(mbathbroadening);
  else
    return -1;
}

LineShape::Output Absorption::Lines::ShapeParameters_dVMR(size_t k, Numeric T, Numeric P, const QuantumIdentifier& vmr_qi) const noexcept {
  const auto self = vmr_qi.Species() == mquantumidentity.Species();
  const auto& ls = mlines[k].LineShape();
  
  if (mselfbroadening and self) {
    auto x = ls.GetVMRDerivs(T, mT0, P, 0);
    
    if (mbathbroadening)
      x = LineShape::differenceOutput(x, ls.GetVMRDerivs(
          T, mT0, P, ls.nelem() - 1));
    
    if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
    
    return x;
  } else if (mbathbroadening and self)
    return {0, 0, 0, 0, 0, 0, 0, 0, 0};
  else {
    auto x = ls.GetVMRDerivs(T, mT0, P, LineShapePos(vmr_qi));
    
    if (mbathbroadening)
      x = LineShape::differenceOutput(x, ls.GetVMRDerivs(
          T, mT0, P, ls.nelem() - 1));
    
    if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
    
    return x;
  }
}

Numeric Absorption::Lines::ShapeParameter_dInternal(size_t k, Numeric T, Numeric P, const Vector& vmrs, const RetrievalQuantity& derivative) const noexcept {
  const auto self = derivative.Mode() == LineShape::self_broadening;
  const auto bath = derivative.Mode() == LineShape::bath_broadening;
  const auto& ls = mlines[k].LineShape();
  
  if(derivative.QuantumIdentity().Species() not_eq Species() or
      derivative.QuantumIdentity().Isotopologue() not_eq Isotopologue())
    return 0;
  else if(self and mselfbroadening)
    return ls.GetInternalDeriv(
      T, mT0, P, 0, vmrs, derivative.PropMatType());
  else if(self)
    return ls.GetInternalDeriv(
      T, mT0, P, LineShapePos(SpeciesTag(derivative.Mode()).Species()), vmrs, derivative.PropMatType());
  else if(bath and mbathbroadening)
    return ls.GetInternalDeriv(
      T, mT0, P, ls.nelem() - 1, vmrs, derivative.PropMatType());
  else if(bath)
    return 0;
  else
    return ls.GetInternalDeriv(
      T, mT0, P, LineShapePos(SpeciesTag(derivative.Mode()).Species()), vmrs, derivative.PropMatType());
}

Absorption::SingleLineInternalRead Absorption::ReadFromArtscat3Stream(istream& is) {
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

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return {SingleLine(), true};

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return {SingleLine(), true};

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

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0) {
    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
    if (i == ArtsMap.end()) {
      ostringstream os;
      os << "ARTS Tag: " << artsid << " is unknown.";
      throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;
    
    QuantumIdentifier quantumidentity;
    quantumidentity.SetTransition();

    // Set mspecies:
    quantumidentity.SetSpecies(id.Speciesindex());

    // Set misotopologue:
    quantumidentity.SetIsotopologue(id.Isotopologueindex());

    // Extract center frequency:
    Numeric F0;
    icecream >> F0;

    Numeric psf;
    // Extract pressure shift:
    icecream >> psf;

    // Extract intensity:
    Numeric I0;
    icecream >> I0;

    // Extract reference temperature for Intensity in K:
    Numeric T0;
    icecream >> T0;

    // Extract lower state energy:
    Numeric E0;
    icecream >> E0;

    // Extract air broadening parameters:
    Numeric agam, sgam;
    icecream >> agam;
    icecream >> sgam;

    // Extract temperature coefficient of broadening parameters:
    Numeric nair, nself;
    icecream >> nair;
    icecream >> nself;

    // Extract reference temperature for broadening parameter in K:
    Numeric tgam;
    icecream >> tgam;

    // Extract the aux parameters:
    Index naux;
    icecream >> naux;

    // resize the aux array and read it
    ArrayOfNumeric maux;
    maux.resize(naux);

    for (Index j = 0; j < naux; j++) {
      icecream >> maux[j];
    }

    // Extract accuracies:
    Numeric dagam, dsgam, dnair, dnself, dpsf;
    try {
      Numeric mdf, mdi0;
      icecream >> mdf;
      icecream >> mdi0;
      icecream >> dagam;
      icecream >> dsgam;
      icecream >> dnair;
      icecream >> dnself;
      icecream >> dpsf;
    } catch (const std::runtime_error&) {
      // Nothing to do here, the accuracies are optional, so we
      // just set them to -1 and continue reading the next line of
      // the catalogue
      dagam = -1;
      dsgam = -1;
      dnair = -1;
      dnself = -1;
      dpsf = -1;
    }

    // Fix if tgam is different from ti0
    if (tgam != T0) {
      agam = agam * pow(tgam / T0, nair);
      sgam = sgam * pow(tgam / T0, nself);
      psf = psf * pow(tgam / T0, (Numeric).25 + (Numeric)1.5 * nair);
    }
    
    return {SingleLine(F0, I0, E0, NAN, NAN, NAN, Zeeman::Model(), LineShape::Model2(sgam, nself, agam, nair, psf)), false};
  }
  throw std::runtime_error("Bad stream");
}

Absorption::SingleLineInternalRead Absorption::ReadFromArtscat4Stream(istream& is, const std::vector<QuantumNumberType>& localquantums) {
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

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return {SingleLine(), true};

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return {SingleLine(), true};

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

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0) {
    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
    if (i == ArtsMap.end()) {
      ostringstream os;
      os << "ARTS Tag: " << artsid << " is unknown.";
      throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;

    // Set line ID
    QuantumIdentifier quantumidentity;
    quantumidentity.SetTransition();
    quantumidentity.SetSpecies(id.Speciesindex());
    quantumidentity.SetIsotopologue(id.Isotopologueindex());

    // Extract center frequency:
    Numeric F0;
    icecream >> F0;

    // Extract intensity:
    Numeric I0;
    icecream >> I0;

    // Extract reference temperature for Intensity in K:
    Numeric T0;
    icecream >> T0;

    // Extract lower state energy:
    Numeric E0;
    icecream >> E0;

    // Extract Einstein A-coefficient:
    Numeric A;
    icecream >> A;

    // Extract upper state stat. weight:
    Numeric gupp;
    icecream >> gupp;

    // Extract lower state stat. weight:
    Numeric glow;
    icecream >> glow;

    LineShape::Model x;
    LineShape::from_artscat4(icecream, x, quantumidentity);

    // Remaining entries are the quantum numbers
    String mquantum_numbers_str;
    getline(icecream, mquantum_numbers_str);

    mquantum_numbers_str.trim();
    // FIXME: OLE: Added this if to catch crash for species like CO, PH3
    // where the line in the catalog is too short. Better would be to
    // only read the n and j for Zeeman species, but we don't have that
    // information here

    if (species_data[quantumidentity.Species()].Name() == "SO") {
      // Note that atoi converts *** to 0.
      quantumidentity.UpperQuantumNumbers().Set(
          QuantumNumberType::N,
          Rational(atoi(mquantum_numbers_str.substr(0, 3).c_str())));
      quantumidentity.LowerQuantumNumbers().Set(
          QuantumNumberType::N,
          Rational(atoi(mquantum_numbers_str.substr(6, 3).c_str())));
      quantumidentity.UpperQuantumNumbers().Set(
          QuantumNumberType::J,
          Rational(atoi(mquantum_numbers_str.substr(3, 3).c_str())));
      quantumidentity.LowerQuantumNumbers().Set(
          QuantumNumberType::J,
          Rational(atoi(mquantum_numbers_str.substr(9, 3).c_str())));
    }

    if (mquantum_numbers_str.nelem() >= 25) {
      if (species_data[quantumidentity.Species()].Name() == "O2") {
        String vstr = mquantum_numbers_str.substr(0, 3);
        ArrayOfIndex v(3);
        for (Index vi = 0; vi < 3; vi++) {
          if (vstr[0] != ' ')
            extract(v[vi], vstr, 1);
          else
            v[vi] = -1;
        }

        if (v[2] > -1) {
          quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::v1, Rational(v[2]));
          quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::v1, Rational(v[2]));
        }

        String qstr1 = mquantum_numbers_str.substr(4, 12);
        String qstr2 = mquantum_numbers_str.substr(4 + 12 + 1, 12);
        ArrayOfIndex q(6);
        for (Index qi = 0; qi < 3; qi++) {
          if (qstr1.substr(0, 4) != "    ")
            extract(q[qi], qstr1, 4);
          else
            q[qi] = -1;
        }
        for (Index qi = 3; qi < 6; qi++) {
          if (qstr2.substr(0, 4) != "    ")
            extract(q[qi], qstr2, 4);
          else
            q[qi] = -1;
        }

        if (q[0] > -1)
          quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::N, Rational(q[0]));
        if (q[1] > -1)
          quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::J, Rational(q[1]));
        if (q[2] > -1)
          quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::F,
                                         q[2] - Rational(1, 2));
        if (q[3] > -1)
          quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::N, Rational(q[3]));
        if (q[4] > -1)
          quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::J, Rational(q[4]));
        if (q[5] > -1)
          quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::F,
                                         q[5] - Rational(1, 2));
      }
    }
    
    std::vector<Rational> lqn(localquantums.size());
    std::vector<Rational> uqn(localquantums.size());
    for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
      lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
      uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
    }
    
    return {SingleLine(F0, I0, E0, glow, gupp, A, Zeeman::Model(), LineShape::Model2(x.Data()), lqn, uqn), false};
  }
  throw std::runtime_error("Bad stream");
}

Absorption::SingleLineInternalRead Absorption::ReadFromArtscat5Stream(istream& is, const std::vector<QuantumNumberType>& localquantums) {
  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  LineShape::Model line_mixing_model;
  bool lmd_found = false;

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

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return {SingleLine(), true};

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return {SingleLine(), true};

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
  
  LineShape::Model x;
  Zeeman::Model zeeman;
  QuantumIdentifier quantumidentity;
  Numeric F0;
  Numeric I0;
  Numeric T0;
  Numeric E0;
  Numeric A;
  Numeric gupp;
  Numeric glow;
  
  try {
    String artsid;
    icecream >> artsid;

    if (artsid.length() != 0) {
      // ok, now for the cool index map:
      // is this arts identifier valid?
      const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
      if (i == ArtsMap.end()) {
        ostringstream os;
        os << "ARTS Tag: " << artsid << " is unknown.";
        throw runtime_error(os.str());
      }

      SpecIsoMap id = i->second;

      // Set line ID:
      quantumidentity.SetTransition();
      quantumidentity.SetSpecies(id.Speciesindex());
      quantumidentity.SetIsotopologue(id.Isotopologueindex());

      // Extract center frequency:
      icecream >> F0;

      // Extract intensity:
      icecream >> I0;

      // Extract reference temperature for Intensity in K:
      icecream >> T0;

      // Extract lower state energy:
      icecream >> E0;

      // Extract Einstein A-coefficient:
      icecream >> A;

      // Extract upper state stat. weight:
      icecream >> gupp;

      // Extract lower state stat. weight:
      icecream >> glow;

      String token;
      Index nelem;
      icecream >> token;

      while (icecream) {
        // Read pressure broadening (LEGACY)
        if (token == "PB") {
          LineShape::from_pressurebroadeningdata(
              icecream, x, quantumidentity);
          icecream >> token;
        } else if (token == "QN") {
          // Quantum numbers

          icecream >> token;
          if (token != "UP") {
            ostringstream os;
            os << "Unknown quantum number tag: " << token;
            throw std::runtime_error(os.str());
          }

          icecream >> token;
          Rational r;
          while (icecream) {
            ThrowIfQuantumNumberNameInvalid(token);
            icecream >> r;
            quantumidentity.UpperQuantumNumbers().Set(token, r);
            icecream >> token;
            if (token == "LO") break;
          }

          if (!is || token != "LO") {
            std::ostringstream os;
            os << "Error in catalog. Lower quantum number tag 'LO' not found.";
            throw std::runtime_error(os.str());
          }

          icecream >> token;
          while (icecream && IsValidQuantumNumberName(token)) {
            icecream >> r;
            quantumidentity.LowerQuantumNumbers().Set(token, r);
            icecream >> token;
          }
        } else if (token == "LM")  // LEGACY
        {
          LineShape::from_linemixingdata(icecream, line_mixing_model);
          icecream >> token;
          lmd_found = true;
        } else if (token == "LF") {
          LineShape::from_linefunctiondata(icecream, x);
          icecream >> token;
        } else if (token == "LS") {
          icecream >> x;
          icecream >> token;
        } else if (token == "ZM") {
          // Zeeman effect
          icecream >> zeeman;
          icecream >> token;
        } else if (token == "LSM") {
          // Line shape modifications

          // Starts with the number of modifications
          icecream >> nelem;
          for (Index lsm = 0; lsm < nelem; lsm++) {
            icecream >> token;

            // cutoff frequency
            if (token == "CUT") {
              Numeric value = NAN;
              icecream >> value;
            }

            // linemixing pressure limit
            if (token == "LML") {
              Numeric value = NAN;
              icecream >> value;
            }

            // mirroring
            else if (token == "MTM") {
              String value;
              icecream >> value;
            }

            // line normalization
            else if (token == "LNT") {
              String value;
              icecream >> value;
            }

            else {
              ostringstream os;
              os << "Unknown line modifications given: " << token;
              throw std::runtime_error(os.str());
            }
          }
          icecream >> token;
        } else {
          ostringstream os;
          os << "Unknown line data tag: " << token;
          throw std::runtime_error(os.str());
        }
      }
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Parse error in catalog line: " << endl;
    os << line << endl;
    os << e.what();
    throw std::runtime_error(os.str());
  }

  if (lmd_found)
    x.SetLineMixingModel(line_mixing_model.Data()[0]);

  std::vector<Rational> lqn(localquantums.size());
  std::vector<Rational> uqn(localquantums.size());
  for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
    lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
    uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
  }
  
  return {SingleLine(F0, I0, E0, glow, gupp, A, Zeeman::Model(), LineShape::Model2(x.Data()), lqn, uqn), false};
}

// The below is a derivative of ReadFromHitran2001Stream
Absorption::SingleLineInternalRead Absorption::ReadFromLBLRTMStream(istream& is, const std::vector<QuantumNumberType>& localquantums) {
  // Global species lookup data:
  using global_data::species_data;

  // This value is used to flag missing data both in species and
  // isotopologue lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotopologue.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>].
  //
  // Allow for up to 100 species in HITRAN in the future.
  static Array<Index> hspec(100);

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopologues.
  static Array<ArrayOfIndex> hiso(100);

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning:
  static ArrayOfIndex warned_missing;

  if (!hinit) {
    // Initialize hspec.
    // The value of missing means that we don't have this species.
    hspec = missing;  // Matpack can set all elements like this.
    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];
      // We have to be careful and check for the case that all
      // HITRAN isotopologue tags are -1 (this species is missing in HITRAN).
      if (sr.Isotopologue().nelem() && 0 < sr.Isotopologue()[0].HitranTag()) {
        // The HITRAN tags are stored as species plus isotopologue tags
        // (MO and ISO)
        // in the Isotopologue() part of the species record.
        // We can extract the MO part from any of the isotopologue tags,
        // so we use the first one. We do this by taking an integer
        // division by 10.

        Index mo = sr.Isotopologue()[0].HitranTag() / 10;
        //          cout << "mo = " << mo << endl;
        hspec[mo] = i;

        // Get a nicer to handle array of HITRAN iso tags:
        Index n_iso = sr.Isotopologue().nelem();
        ArrayOfIndex iso_tags;
        iso_tags.resize(n_iso);
        for (Index j = 0; j < n_iso; ++j) {
          iso_tags[j] = sr.Isotopologue()[j].HitranTag();
        }

        // Reserve elements for the isotopologue tags. How much do we
        // need? This depends on the largest HITRAN tag that we know
        // about!
        // Also initialize the tags to missing.
        //          cout << "iso_tags = " << iso_tags << endl;
        //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
        //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
        hiso[mo].resize(max(iso_tags) % 10 + 1);
        hiso[mo] = missing;  // Matpack can set all elements like this.

        // Set the isotopologue tags:
        for (Index j = 0; j < n_iso; ++j) {
          if (0 < iso_tags[j])  // ignore -1 elements
          {
            // To get the iso tags from HitranTag() we also have to take
            // modulo 10 to get rid of mo.
            hiso[mo][iso_tags[j] % 10] = j;
          }
        }
      }
    }

    hinit = true;
  }

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return {SingleLine(), true};

    // Throw runtime_error if stream is bad:
    if (!is) throw std::runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return {SingleLine(), true};

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    // Because of the fixed FORTRAN format, we need to break up the line
    // explicitly in apropriate pieces. Not elegant, but works!

    // Extract molecule number:
    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    //      cout << "mo = " << mo << endl;

    // If mo == 0 this is just a comment line:
    if (0 != mo) {
      // See if we know this species. Exit with an error if the species is unknown.
      if (missing != hspec[mo]) {
        comment = false;

        // Check if data record has the right number of characters for the
        // in Hitran 1986-2001 format
        Index nChar = line.nelem() + 2;  // number of characters in data record;
        if (nChar != 100) {
          ostringstream os;
          os << "Invalid HITRAN 1986-2001 line data record with " << nChar
             << " characters (expected: 100)." << endl
             << line << " n: " << line.nelem();
          throw std::runtime_error(os.str());
        }

      }
    }
  }
  
  LineShape::Model x;
  QuantumIdentifier quantumidentity;
  Numeric F0;
  Numeric I0;
  Numeric E0;

  // Set mspecies from my cool index table:
  quantumidentity.SetTransition();
  quantumidentity.SetSpecies(hspec[mo]);

  // Extract isotopologue:
  Index iso;
  extract(iso, line, 1);
  //  cout << "iso = " << iso << endl;

  // Set misotopologue from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases.
  quantumidentity.SetIsotopologue(missing);
  if (iso < hiso[mo].nelem())
    if (missing != hiso[mo][iso]) quantumidentity.SetIsotopologue(hiso[mo][iso]);

  // Issue error message if misotopologue is still missing:
  if (missing == quantumidentity.Isotopologue()) {
    ostringstream os;
    os << "Species: " << species_data[quantumidentity.Species()].Name()
       << ", isotopologue iso = " << iso << " is unknown.";
    throw std::runtime_error(os.str());
  }

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    F0 = v * w2Hz;
    //    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT;  // in [m/s]

    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;
    if (line[6] == 'D') line[6] = 'E';
    // Extract HITRAN intensity:
    extract(s,
            line,
            10);  // NOTE:  If error shooting, FORTRAN "D" is not read properly.
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    I0 = s * hi2arts;
    // Take out isotopologue ratio:
    I0 /= species_data[quantumidentity.Species()]
               .Isotopologue()[quantumidentity.Isotopologue()]
               .Abundance();
  }

  // Skip transition probability:
  {
    Numeric r;
    extract(r, line, 10);
  }

  // Air broadening parameters.
  Numeric sgam, agam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(E0, line, 10);

    // Convert to Joule:
    E0 = wavenumber_to_joule(E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  //Skip upper state global quanta index
  {
    Index eu;
    extract(eu, line, 3);
  }

  //Skip lower state global quanta index
  {
    Index el;
    extract(el, line, 3);
  }

  //Skip upper state local quanta
  {
    Index eul;
    extract(eul, line, 9);
  }

  //Skip lower state local quanta
  {
    Index ell;
    if (species_data[quantumidentity.Species()].Name() == "O2") {
      String helper = line.substr(0, 9);
      Index DJ = -helper.compare(3, 1, "Q");
      Index DN = -helper.compare(0, 1, "Q");
      Index N = atoi(helper.substr(1, 2).c_str());
      Index J = atoi(helper.substr(4, 2).c_str());

      quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::N, N);
      quantumidentity.LowerQuantumNumbers().Set(QuantumNumberType::J, J);
      quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::N, N - DN);
      quantumidentity.UpperQuantumNumbers().Set(QuantumNumberType::J, J - DJ);
    }

    extract(ell, line, 9);
  }

  // Accuracy index for frequency reference
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity reference
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for halfwidth reference
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Skip four
  {
    Index four;
    extract(four, line, 4);
  }

  // This is the test for the last two characters of the line
  {
    /* 
     *   0 is nothing, 
     *  -1 is linemixing on the next line, 
     *  -3 is the non-resonant line 
     */
    Index test;
    extract(test, line, 2);
    //If the tag is as it should be, then a minus one means that more should be read
    if (test == -1 || test == -3)
      getline(is, line);
    else  // the line is done and we are happy to leave
    {
      x = LineShape::Model(sgam, nself, agam, nair, psf);
      
      std::vector<Rational> lqn(localquantums.size());
      std::vector<Rational> uqn(localquantums.size());
      for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
        lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
        uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
      }
      
      return {SingleLine(F0, I0, E0, NAN, NAN, NAN, Zeeman::Model(), LineShape::Model2(x.Data()), lqn, uqn), false};
    }
  }

  // In case we are unable to leave, the next line is a line mixing parameter line

  // First is the molecular number.  This should be the same as above.
  {
    Index mo2;
    extract(mo2, line, 2);
    // Skip one

    if (mo != mo2)
      throw std::runtime_error("There is an error in the line mixing\n");
  }

  Vector Y(4), G(4), T(4);

  // These are constants for AER but should be included because we need their grid.
  T[0] = 200;
  T[1] = 250;
  T[2] = 296;
  T[3] = 340;

  // Next is the Y  and G at various temperatures
  {
    Numeric Y_200K;
    extract(Y_200K, line, 13);
    Y[0] = Y_200K;
  }
  {
    Numeric G_200K;
    extract(G_200K, line, 11);
    G[0] = G_200K;
  }
  {
    Numeric Y_250K;
    extract(Y_250K, line, 13);
    Y[1] = Y_250K;
  }
  {
    Numeric G_250K;
    extract(G_250K, line, 11);
    G[1] = G_250K;
  }
  {
    Numeric Y_296K;
    extract(Y_296K, line, 13);
    Y[2] = Y_296K;
  }
  {
    Numeric G_296K;
    extract(G_296K, line, 11);
    G[2] = G_296K;
  }
  {
    Numeric Y_340K;
    extract(Y_340K, line, 13);
    Y[3] = Y_340K;
  }
  {
    Numeric G_340K;
    extract(G_340K, line, 11);
    G[3] = G_340K;
  }

  extern const Numeric ATM2PA;

  Y /= ATM2PA;
  G /= ATM2PA / ATM2PA;
  Y *=
      -1;  // ARTS uses (1-iY) as line-mixing factor, LBLRTM CO2 uses (1+iY), so we must change sign

  // Test that this is the end
  {
    Index test;
    extract(test, line, 2);
    if (test == -1) {
      x = LineShape::Model(sgam,
                           nself,
                           agam,
                           nair,
                           psf,
                           {T[0],
                            T[1],
                            T[2],
                            T[3],
                            Y[0],
                            Y[1],
                            Y[2],
                            Y[3],
                            G[0],
                            G[1],
                            G[2],
                            G[3]});
      
      std::vector<Rational> lqn(localquantums.size());
      std::vector<Rational> uqn(localquantums.size());
      for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
        lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
        uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
      }
      
      return {SingleLine(F0, I0, E0, NAN, NAN, NAN, Zeeman::Model(), LineShape::Model2(x.Data()), lqn, uqn), false};
    } else if (test == -3) {
      x = LineShape::Model(sgam,
                           nself,
                           agam,
                           nair,
                           psf,
                           {T[0],
                            T[1],
                            T[2],
                            T[3],
                            Y[0],
                            Y[1],
                            Y[2],
                            Y[3],
                            G[0],
                            G[1],
                            G[2],
                            G[3]});
      
      std::vector<Rational> lqn(localquantums.size());
      std::vector<Rational> uqn(localquantums.size());
      for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
        lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
        uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
      }
      
      return {SingleLine(F0, I0, E0, NAN, NAN, NAN, Zeeman::Model(), LineShape::Model2(x.Data()), lqn, uqn), false};
    } else
      return {SingleLine(), true};
  }
}

Absorption::SingleLineInternalRead Absorption::ReadFromHitran2004Stream(istream& is, const std::vector<QuantumNumberType>& localquantums, Numeric fmin) {
  // Global species lookup data:
  using global_data::species_data;

  // This value is used to flag missing data both in species and
  // isotopologue lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotopologue.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>].
  //
  // Allow for up to 100 species in HITRAN in the future.
  static Array<Index> hspec(100);

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopologues.
  static Array<ArrayOfIndex> hiso(100);

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning:
  static ArrayOfIndex warned_missing;

  if (!hinit) {
    // Initialize hspec.
    // The value of missing means that we don't have this species.
    hspec = missing;  // Matpack can set all elements like this.

    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];

      // We have to be careful and check for the case that all
      // HITRAN isotopologue tags are -1 (this species is missing in HITRAN).

      if (sr.Isotopologue().nelem() && 0 < sr.Isotopologue()[0].HitranTag()) {
        // The HITRAN tags are stored as species plus isotopologue tags
        // (MO and ISO)
        // in the Isotopologue() part of the species record.
        // We can extract the MO part from any of the isotopologue tags,
        // so we use the first one. We do this by taking an integer
        // division by 10.

        Index mo = sr.Isotopologue()[0].HitranTag() / 10;
        //          cout << "mo = " << mo << endl;
        hspec[mo] = i;

        // Get a nicer to handle array of HITRAN iso tags:
        Index n_iso = sr.Isotopologue().nelem();
        ArrayOfIndex iso_tags;
        iso_tags.resize(n_iso);
        for (Index j = 0; j < n_iso; ++j) {
          iso_tags[j] = sr.Isotopologue()[j].HitranTag();
        }

        // Reserve elements for the isotopologue tags. How much do we
        // need? This depends on the largest HITRAN tag that we know
        // about!
        // Also initialize the tags to missing.
        //          cout << "iso_tags = " << iso_tags << endl;
        //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
        //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
        hiso[mo].resize(max(iso_tags) % 10 + 1);
        hiso[mo] = missing;  // Matpack can set all elements like this.

        // Set the isotopologue tags:
        for (Index j = 0; j < n_iso; ++j) {
          if (0 < iso_tags[j])  // ignore -1 elements
          {
            // To get the iso tags from HitranTag() we also have to take
            // modulo 10 to get rid of mo.
            hiso[mo][iso_tags[j] % 10] = j;
          }
        }
      }
    }

    hinit = true;
  }

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return {SingleLine(), true};

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return {SingleLine(), true};

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    // Because of the fixed FORTRAN format, we need to break up the line
    // explicitly in apropriate pieces. Not elegant, but works!

    // Extract molecule number:
    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    //      cout << "mo = " << mo << endl;

    // If mo == 0 this is just a comment line:
    if (0 != mo) {
      // See if we know this species.
      if (missing != hspec[mo]) {
        comment = false;

        // Check if data record has the right number of characters for the
        // in Hitran 2004 format
        Index nChar = line.nelem() + 2;  // number of characters in data record;
        if ((nChar == 161 && line[158] != ' ') || nChar > 161) {
          ostringstream os;
          os << "Invalid HITRAN 2004 line data record with " << nChar
             << " characters (expected: 160).";
          throw std::runtime_error(os.str());
        }

      }
    }
  }

  // Ok, we seem to have a valid species here.
  
  LineShape::Model x;
  QuantumIdentifier quantumidentity;
  Numeric F0;
  Numeric I0;
  Numeric E0;
  Numeric A;
  Numeric gupp;
  Numeric glow;

  // Set mspecies from my cool index table:
  quantumidentity.SetTransition();
  quantumidentity.SetSpecies(hspec[mo]);

  // Extract isotopologue:
  Index iso;
  extract(iso, line, 1);
  //  cout << "iso = " << iso << endl;

  // Set misotopologue from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases.
  quantumidentity.SetIsotopologue(missing);
  if (iso < hiso[mo].nelem())
    if (missing != hiso[mo][iso]) quantumidentity.SetIsotopologue(hiso[mo][iso]);

  // Issue error message if misotopologue is still missing:
    if (missing == quantumidentity.Isotopologue()) {
    ostringstream os;
    os << "Species: " << species_data[quantumidentity.Species()].Name()
       << ", isotopologue iso = " << iso << " is unknown.";
    throw std::runtime_error(os.str());
  }

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    F0 = v * w2Hz;
    //    cout << "mf = " << mf << endl;
    if (F0 < fmin) {
      F0 = -1;
      return {SingleLine(), false};
    }
  }

  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT;  // in [m/s]

    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s, line, 10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    I0 = s * hi2arts;
    // Take out isotopologue ratio:
    I0 /= species_data[quantumidentity.Species()]
              .Isotopologue()[quantumidentity.Isotopologue()]
              .Abundance();
  }

  // Einstein coefficient
  {
    Numeric r;
    extract(r, line, 10);
    A = r;
  }

  // Air broadening parameters.
  Numeric agam, sgam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(E0, line, 10);

    // Convert to Joule:
    E0 = wavenumber_to_joule(E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  static QuantumParserHITRAN2004 quantum_parser;
  const String qstr = line.substr(0, 15 * 4);

  // Upper state global quanta
  {
    Index eu;
    extract(eu, line, 15);
  }

  // Lower state global quanta
  {
    Index el;
    extract(el, line, 15);
  }

  // Upper state local quanta
  {
    Index eul;
    extract(eul, line, 15);
  }

  // Lower state local quanta
  {
    Index ell;
    extract(ell, line, 15);
  }

  // Parse quantum numbers.
  quantum_parser.Parse(quantumidentity, qstr);

  // Accuracy index for frequency
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for air-broadened halfwidth
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for self-broadened half-width
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for temperature-dependence exponent for agam
  {
    Index dn;
    // Extract HITRAN value:
    extract(dn, line, 1);
  }

  // Accuracy index for pressure shift
  {
    Index dpsfi;
    // Extract HITRAN value (given in cm-1):
    extract(dpsfi, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN 2004. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Set line shape computer
  x = LineShape::Model(sgam, nself, agam, nair, psf);
  {
    Index garbage;
    extract(garbage, line, 13);

    // The statistical weights
    extract(gupp, line, 7);
    extract(glow, line, 7);
  }
  
  std::vector<Rational> lqn(localquantums.size());
  std::vector<Rational> uqn(localquantums.size());
  for(size_t kqn=0; kqn<localquantums.size(); kqn++) {
    lqn[kqn] = quantumidentity.LowerQuantumNumber(QuantumNumberType(kqn));
    uqn[kqn] = quantumidentity.UpperQuantumNumber(QuantumNumberType(kqn));
  }
  
  // That's it!
  return {SingleLine(F0, I0, E0, glow, gupp, A, Zeeman::Model(), LineShape::Model2(sgam, nself, agam, nair, psf), lqn, uqn), false};
}
