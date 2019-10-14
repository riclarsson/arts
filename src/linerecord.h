/* Copyright (C) 2000-2013
   Stefan Buehler <sbuehler@ltu.se>
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

/** \file
    LineRecord class for managing line catalog data.

    \author Stefan Buehler, Axel von Engeln
*/

#ifndef linerecord_h
#define linerecord_h

#include <cmath>
#include <stdexcept>
#include "array.h"
#include "constants.h"
#include "lineshapemodel.h"
#include "matpackI.h"
#include "messages.h"
#include "mystring.h"
#include "quantum.h"
#include "wigner_functions.h"
#include "zeemandata.h"

/* Forward declaration of classes */
class SpeciesRecord;
class IsotopologueRecord;

/** Spectral line catalog data.

    Below is a description of the ARTS catalogue format.
    The file starts with the usual XML header:

    \verbatim
    <?xml version="1.0"?>
    <arts format="ascii" version="1">
    <ArrayOfLineRecord version="ARTSCAT-3" nelem="8073">
    \endverbatim

    The ARTSCAT version number is there to keep track of catalogue format changes.
    The "nelem" tag contains the total number of lines in the file.

    The file ends with the usual XML closing tags:

    \verbatim
    </ArrayOfLineRecord>
    </arts>
    \endverbatim

    In-between the header and the footer are the actual spectroscopic
    data. Each new entry, corresponding to one spectral line, starts
    with the `@' character.

    The line catalogue should not have any fixed column widths because
    the precision of the parameters should not be limited by the
    format.  The catalogue can then be stored as binary or ASCII. In
    the ASCII version the columns are separated by one or more
    blanks. The line format is then specified by only the order and
    the units of the columns. As the catalogue entry for each
    transition can be quite long, it can be broken across lines in the
    ASCII file. That is why each new transition is marked with a `@'
    character.

    The first column will contain the species and isotopologue, following
    the naming scheme described below.  Scientific notation is
    allowed, e.g. 501.12345e9.  

    Note that starting with ARTSCAT-2, the intensity is per molecule,
    i.e., it does not contain the isotopologue ratio. This is similar to
    JPL, but different to HITRAN.

    Currently, ARTS is capable of handling ARTSCAT versions 3 and 4. Different
    versions can be handled simultaneously.
    The line format of ARTSCAT-3 (for ARTSCAT-4 see further below) is:

    \verbatim
    Col  Variable                    Label     Unit    Comment
    ------------------------------------------------------------------      
     0   `@'                         ENTRY        -    marks start of entry
     1    species\&isotopologue tag   NAME        -    e.g. O3-666
     2   center frequency                F       Hz    e.g. 501.12345e9 
     3   pressure shift of F           PSF    Hz/Pa    
     4   line intensity                 I0    Hz*m^2   per isotopologue, not per species
     5   reference temp. for I0       T_I0        K
     6   lower state energy           ELOW        J    
     7   air broadened width          AGAM    Hz/Pa    values around 20000 Hz/Pa
     8   self broadened width         SGAM    Hz/Pa
     9   AGAM temp. exponent          NAIR        -    values around .5
    10   SGAM temp. exponent         NSELF        - 
    11   ref. temp. for AGAM, SGAM   T_GAM        K
    12   number of aux. parameters   N_AUX        -
    13   auxiliary parameter          AUX1        -
    14   ... 
    15   error for F                    DF       Hz
    16   error for I0                  DI0        %
    17   error for AGAM              DAGAM        %
    18   error for SGAM              DSGAM        %
    19   error for NAIR              DNAIR        %
    20   error for NSELF            DNSELF        %
    21   error for PSF                DPSF        %
    \endverbatim

    The parameters 0-12 must be present, the others can be missing,
    since they are not needed for the calculation.

    For the error fields (15-21), a -1 means that no value exist.

    A valid ARTS (CAT-3) line file would be:
    \verbatim
    <?xml version="1.0"?>
    <arts format="ascii" version="1">
    <ArrayOfLineRecord version="ARTSCAT-3" nelem="2">
    @ O3-676 80015326542.0992 0 3.70209114155527e-19 296 7.73661776567701e-21 21480.3182341969
    28906.7092490501 0.76 0.76 296 0 300000 0.1 0.1 0.1 0.2 -1 -1.24976056865038e-11
    @ O3-676 80015476438.3282 0 3.83245786810611e-19 296 7.73661577919822e-21 21480.3182341969
    28906.7092490501 0.76 0.76 296 0 300000 0.1 0.1 0.1 0.2 -1 -1.24975822742335e-11
    </ArrayOfLineRecord>
    </arts>
    \endverbatim

    Some species need special parameters that are not needed by other
    species (for example overlap coefficients for O2). In the case of
    oxygen two parameters are sufficient to describe the overlap, but
    other species, e.g., methane, may need more coefficients. The
    default for \c N_AUX is zero. In that case, no further
    \c AUX fields are present. [FIXME: Check Oxygen.]

    The line format of ARTSCAT-4 is:

    \verbatim
    Col  Variable                    Label        Unit    Comment
    ------------------------------------------------------------------      
    00   `@'                         ENTRY        -       marks start of entry
    01   species\&isotopologue tag   NAME         -       e.g. O3-666
    02   center frequency            F            Hz      e.g. 501.12345e9 
    03   line intensity              I0           Hz*m^2  per isotopologue, not per species
    04   reference temp. for I0      T_I0         K
    05   lower state energy          ELOW         J   
    06   Einstein A-coefficient      A            1/s     where available from HITRAN 
    07   Upper state stat. weight    G_upper      -       where available from HITRAN 
    08   Lower state stat. weight    G_lower      -       where available from HITRAN 

    09   broadening parameter self   GAMMA_self   Hz/Pa   GAM have values around 20000 Hz/Pa
    10   broadening parameter N2     GAMMA_N2     Hz/Pa
    11   broadening parameter O2     GAMMA_O2     Hz/Pa
    12   broadening parameter H2O    GAMMA_H2O    Hz/Pa
    13   broadening parameter CO2    GAMMA_CO2    Hz/Pa
    14   broadening parameter H2     GAMMA_H2     Hz/Pa
    15   broadening parameter He     GAMMA_He     Hz/Pa

    16   GAM temp. exponent self     N_self       -       N have values around .5
    17   GAM temp. exponent N2       N_N2         -       
    18   GAM temp. exponent O2       N_O2         -       
    19   GAM temp. exponent H2O      N_H2O        -       
    20   GAM temp. exponent CO2      N_CO2        -       
    21   GAM temp. exponent H2       N_H2         -       
    22   GAM temp. exponent He       N_He         -       

    23   F pressure shift N2         DELTA_N2     Hz/Pa   DELTA have values around 0 Hz/Pa
    24   F pressure shift O2         DELTA_O2     Hz/Pa    
    25   F pressure shift H2O        DELTA_H2O    Hz/Pa    
    26   F pressure shift CO2        DELTA_CO2    Hz/Pa    
    27   F pressure shift H2         DELTA_H2     Hz/Pa    
    28   F pressure shift He         DELTA_He     Hz/Pa    

    29   Vibrational and rotational  VRA          -       contains (coded) quantum numbers.
         assignments                                     
------------------------------------------------------------------      
\endverbatim

    Parameters 0-28 must be present.
    Coding conventions of parameter 29 are species specific. The definition is
    given in arts-xml-data/spectroscopy/perrin/ARTSCAT-4_Col29_Conventions.txt.


    The names of the private members and public access functions of
    this data structure follow the above table. The only difference is
    that underscores are omited and only the first letter of each name
    is capitalized. This is for consistency with the notation
    elsewhere in the program.

    \author Stefan Buehler 
*/

enum class MirroringType : Index {
  // Describes the type of mirroring line effects.
  // LineRecord effect:  Adds MTM (String)MirroringType under line modifications when saving ARTSCAT5
  None,             // No mirroring
  Lorentz,          // Mirror, but use Lorentz line shape
  SameAsLineShape,  // Mirror using the same line shape
  Manual,  // Mirror by having a line in the array of line record with negative F0
};

enum class LineNormalizationType : Index {
  // Describes the type of normalization line effects.
  // LineRecord effect:  Adds LNT (String)LineNormalizationType under line modifications when saving ARTSCAT5
  // 0 is ignored on printing
  None,  // Do not renormalize the line shape
  VVH,   // Renormalize with Van Vleck and Huber specifications
  VVW,   // Renormalize with Van Vleck and Weiskopf specifications
  RosenkranzQuadratic,  // Renormalize using Rosenkranz's quadratic specifications
};

enum class LinePopulationType : Index {
  // Describes the type of population level counter.
  // LineRecord effect:  Adds POP (String)LinePopulationType under line modifications when saving ARTSCAT5
  ByLTE,                      // Assume line is in LTE
  ByVibrationalTemperatures,  // Assume line is in NLTE described by vibrational temperatures
  ByPopulationDistribution,  // Assume line is in NLTE and the upper-to-lower ratio is known
};

class LineRecord {
 public:
  /** Default constructor. Initialize to default values. The indices
      are initialized to large numbers, so that we at least get range
      errors when we try to used un-initialized data. */
  LineRecord()
      : mversion(5),
        mqid(QuantumIdentifier::TRANSITION, 1000000, 1000000),
        mf(0.),
        mi0(0.),
        mti0(0.),
        melow(0.),
        mevlow(-1.0),
        mnlte_lower_index(-1),
        mevupp(-1.0),
        mnlte_upper_index(-1),
        ma(NAN),
        mgupper(NAN),
        mglower(NAN),
        mstandard(true),
        mcutoff(-1.0),
        mlinemixing_limit(-1.0),
        mmirroring(MirroringType::None),
        mlinenorm(LineNormalizationType::None),
        mpopulation(LinePopulationType::ByLTE) { /* Nothing to do here. */
  }

  /** Constructor that sets all data elements explicitly. If
      assertions are not disabled (i.e., if NDEBUG is not \#defined),
      assert statements check that the species and isotopologue data
      exists. */
  LineRecord(Index species,
             Index isotopologue,
             Numeric f,
             Numeric psf,
             Numeric i0,
             Numeric ti0,
             Numeric elow,
             Numeric agam,
             Numeric sgam,
             Numeric nair,
             Numeric nself,
             Numeric /*tgam*/,
             const ArrayOfNumeric& /*aux*/,
             Numeric /* df */,
             Numeric /* di0 */,
             Numeric /* dagam */,
             Numeric /* dsgam */,
             Numeric /* dnair */,
             Numeric /* dnself */,
             Numeric /* dpsf */)
      : mversion(5),
        mqid(QuantumIdentifier::TRANSITION, species, isotopologue),
        mf(f),
        mi0(i0),
        mti0(ti0),
        melow(elow),
        mevlow(-1.0),
        mnlte_lower_index(-1),
        mevupp(-1.0),
        mnlte_upper_index(-1),
        ma(NAN),
        mgupper(NAN),
        mglower(NAN),
        mlineshapemodel(sgam, nself, agam, nair, psf),
        mstandard(true),
        mcutoff(-1.0),
        mlinemixing_limit(-1.0),
        mmirroring(MirroringType::None),
        mlinenorm(LineNormalizationType::None),
        mpopulation(LinePopulationType::ByLTE) {}

  /** Return the version String. */
  String VersionString() const;

  /** Return the version number. */
  Index Version() const { return mversion; }
  void SetVersion5() { mversion = 5; }

  /* Set to latest compatible version of the catalog to store all data */
  void SetVersionToLatest() { mversion = 5; }

  /** The index of the molecular species that this line belongs to.
   The species data can be accessed by species_data[Species()]. */
  Index Species() const { return mqid.Species(); }

  /** The index of the isotopologue species that this line belongs to.
   The isotopologue species data can be accessed by
   species_data[Species()].Isotopologue()[Isotopologue()].  */
  Index Isotopologue() const { return mqid.Isotopologue(); }

  /** The full name of the species and isotopologue. E.g., `O3-666'.
   The name is found by looking up the information in species_data,
   using the species and isotopologue index. */
  String Name() const;

  /** The matching SpeciesRecord from species_data. To get at the
   species data of a LineRecord lr, you can use:
   <ul>
   <li>species_data[lr.Species()]</li>
   <li>lr.SpeciesData()</li>
   </ul>
   The only advantages of the latter are that the notation is
   slightly nicer and that you don't have to declare the external
   variable species_data. */
  const SpeciesRecord& SpeciesData() const;

  /** The matching IsotopologueRecord from species_data. The IsotopologueRecord
   is a subset of the SpeciesRecord. To get at the isotopologue data of
   a LineRecord lr, you can use:
   <ul>
   <li>species_data[lr.Species()].Isotopologue()[lr.Isotopologue()]</li>
   <li>lr.SpeciesData().Isotopologue()[lr.Isotopologue()]</li>
   <li>lr.IsotopologueData()</li>
   </ul>
   The last option is clearly the shortest, and has the advantage
   that you don't have to declare the external variable
   species_data. */
  const IsotopologueRecord& IsotopologueData() const;

  /** The line center frequency in <b> Hz</b>. */
  Numeric F() const { return mf; }

  /** Set the line center frequency in <b> Hz</b>. */
  void setF(Numeric new_mf) { mf = new_mf; }

  /** The pressure shift parameter in <b> Hz/Pa</b>. */
  Numeric Psf() const { return mlineshapemodel.Data().back().D0().X0; }

  /** The line intensity in <b> m^2*Hz</b> at the reference temperature \c Ti0. 

   The line intensity \f$I_0\f$ is defined by:

   \f[
   \alpha(\nu) = n \, x \, I_0(T) \, F(\nu)
   \f]

   where \f$\alpha\f$ is the absorption coefficient (in <b>
   m^-1</b>), \f$\nu\f$ is frequency, \f$n\f$ is the
   total number density, \f$x\f$ is the volume mixing ratio, and
   \f$F(\nu)\f$ is the lineshape function. */
  Numeric I0() const { return mi0; }

  /** Set Intensity */
  void setI0(Numeric new_mi0) { mi0 = new_mi0; }

  /** Reference temperature for I0 in <b> K</b>: */
  Numeric Ti0() const { return mti0; }

  /** Lower state energy in <b> cm^-1</b>: */  //FIXME: really in cm-1 still?
  Numeric Elow() const { return melow; }
  void SetElow(Numeric new_melow) { melow = new_melow; }

  /** Lower state vibrational energy in <b> cm^-1</b>: */  //FIXME: really in cm-1 still?
  Numeric Evlow() const { return mevlow; }
  void SetEvlow(Numeric evlow) { mevlow = evlow; }
  Index NLTELowerIndex() const { return mnlte_lower_index; }
  void SetNLTELowerIndex(Index nlte_lower_index) {
    mnlte_lower_index = nlte_lower_index;
  }

  /** Upper state vibrational energy in <b> cm^-1</b>: */  //FIXME: really in cm-1 still?
  Numeric Evupp() const { return mevupp; }
  void SetEvupp(Numeric evupp) { mevupp = evupp; }
  Index NLTEUpperIndex() const { return mnlte_upper_index; }
  void SetNLTEUpperIndex(Index nlte_upper_index) {
    mnlte_upper_index = nlte_upper_index;
  }

  /** Air broadened width in <b> Hz/Pa</b>: */
  Numeric Agam() const { return mlineshapemodel.Data().back().G0().X0; }

  /** Self broadened width in <b> Hz/Pa</b>: */
  Numeric Sgam() const { return mlineshapemodel.Data().front().G0().X0; }

  /** AGAM temperature exponent (dimensionless): */
  Numeric Nair() const { return mlineshapemodel.Data().back().G0().X1; }

  /** SGAM temperature exponent (dimensionless): */
  Numeric Nself() const { return mlineshapemodel.Data().front().G0().X1; }

  /** ARTSCAT-4/5 Einstein A-coefficient in <b> 1/s </b>: */
  Numeric A() const { return ma; }

  Numeric electric_dipole_moment_squared() const {
    using namespace Constant;
    return 3 * h * epsilon_0 / (16 * pi) * ma * pow3(c / mf);
  }

  Numeric magnetic_quadrapole_moment_squared() const {
    using namespace Constant;
    return 3 * h / (16 * pi * mu_0) * ma * pow3(c / mf);
  }

  Numeric reduced_rovibrational_dipole(Rational k = 1) const {
    const Rational Jf = LowerQuantumNumber(QuantumNumberType::J);
    const Rational Ji = UpperQuantumNumber(QuantumNumberType::J);
    const Rational lf = LowerQuantumNumber(QuantumNumberType::l2);
    const Rational li = UpperQuantumNumber(QuantumNumberType::l2);
    const Numeric val =
        sqrt(2 * Jf + 1) * wigner3j(Jf, k, Ji, li, lf - li, -lf);
    if ((Jf + lf + 1) % 2)
      return -val;
    else
      return +val;
  }

  Numeric reduced_magnetic_quadrapole() const {
    const Rational Jf = LowerQuantumNumber(QuantumNumberType::J);
    const Rational Ji = UpperQuantumNumber(QuantumNumberType::J);
    const Rational N = LowerQuantumNumber(QuantumNumberType::N);
    constexpr Rational one(1, 1);
    const Numeric val = sqrt(6 * (2 * Jf + 1) * (2 * Ji + 1)) *
                        wigner6j(one, one, one, Ji, Jf, N);
    if ((Jf + N) % 2)
      return -val;
    else
      return +val;
  }

  /** ARTSCAT-4/5 Upper state stat. weight: */
  Numeric G_upper() const { return mgupper; }

  /** ARTSCAT-4/5 Lower state stat. weight: */
  Numeric G_lower() const { return mglower; }

  /** Quantum numbers */
  void SetQuantumNumberLower(const Index i, const Rational r) {
    mqid.LowerQuantumNumbers().Set(i, r);
  }
  void SetQuantumNumberLower(const String i, const Rational r) {
    mqid.LowerQuantumNumbers().Set(i, r);
  }
  void SetQuantumNumberLower(const QuantumNumberType i, const Rational r) {
    mqid.LowerQuantumNumbers().Set(i, r);
  }
  void SetQuantumNumberUpper(const QuantumNumberType i, const Rational r) {
    mqid.UpperQuantumNumbers().Set(i, r);
  }
  void SetQuantumNumberUpper(const String i, const Rational r) {
    mqid.UpperQuantumNumbers().Set(i, r);
  }
  void SetQuantumNumberUpper(const Index i, const Rational r) {
    mqid.UpperQuantumNumbers().Set(i, r);
  }

  /** Quantum identifier */
  const QuantumIdentifier& QuantumIdentity() const { return mqid; }
  const QuantumNumbers& LowerQuantumNumbers() const {
    return mqid.LowerQuantumNumbers();
  }
  const QuantumNumbers& UpperQuantumNumbers() const {
    return mqid.UpperQuantumNumbers();
  }
  Rational LowerQuantumNumber(QuantumNumberType X) const {
    return mqid.LowerQuantumNumber(X);
  }
  Rational UpperQuantumNumber(QuantumNumberType X) const {
    return mqid.UpperQuantumNumber(X);
  }
  bool InQuantumID(const QuantumIdentifier& qid) const { return mqid.In(qid); }
  bool UpperStateInQuantumID(const QuantumIdentifier& qid) const {
    return mqid.InUpper(qid);
  }
  bool LowerStateInQuantumID(const QuantumIdentifier& qid) const {
    return mqid.InLower(qid);
  }

  /**Do linemixing test*/
  bool do_linemixing(const Numeric& P) const noexcept {
    if (mlinemixing_limit < 0)
      return true;
    else
      return mlinemixing_limit > P;
  }

  void SetLineMixingParameters(
      Numeric& Y,
      Numeric& G,
      Numeric& DV,
      const Numeric& T,
      const Numeric& P,
      const ConstVectorView vmrs,
      const ArrayOfArrayOfSpeciesTag& abs_species) const {
    if (do_linemixing(P)) {
      const auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
      Y = mlineshapemodel.Y(T, mti0, P, v);
      G = mlineshapemodel.G(T, mti0, P, v);
      DV = mlineshapemodel.DV(T, mti0, P, v);
    } else
      Y = G = DV = 0;
  }

  void SetPressureBroadeningParameters(
      Numeric& G0,
      Numeric& G2,
      Numeric& ETA,
      Numeric& D0,
      Numeric& D2,
      Numeric& FVC,
      const Numeric& T,
      const Numeric& P,
      const ConstVectorView vmrs,
      const ArrayOfArrayOfSpeciesTag& abs_species) const {
    const auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
    G0 = mlineshapemodel.G0(T, mti0, P, v);
    D0 = mlineshapemodel.D0(T, mti0, P, v);
    G2 = mlineshapemodel.G2(T, mti0, P, v);
    D2 = mlineshapemodel.D2(T, mti0, P, v);
    FVC = mlineshapemodel.FVC(T, mti0, P, v);
    ETA = mlineshapemodel.ETA(T, mti0, P, v);
  }

  void SetPressureBroadeningParametersTemperatureDerivative(
      Numeric& dG0,
      Numeric& dG2,
      Numeric& dETA,
      Numeric& dD0,
      Numeric& dD2,
      Numeric& dFVC,
      const Numeric& T,
      const Numeric& P,
      const ConstVectorView vmrs,
      const ArrayOfArrayOfSpeciesTag& abs_species) const {
    const auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
    dG0 = mlineshapemodel.dG0_dT(T, mti0, P, v);
    dD0 = mlineshapemodel.dD0_dT(T, mti0, P, v);
    dG2 = mlineshapemodel.dG2_dT(T, mti0, P, v);
    dD2 = mlineshapemodel.dD2_dT(T, mti0, P, v);
    dFVC = mlineshapemodel.dFVC_dT(T, mti0, P, v);
    dETA = mlineshapemodel.dETA_dT(T, mti0, P, v);
  }

  void SetLineMixing2SecondOrderData(ConstVectorView d) {
    mlineshapemodel.SetLineMixingModel(
        LineShape::LegacyLineMixingData::vector2modellm(
            d, LineShape::LegacyLineMixingData::TypeLM::LM_2NDORDER)
            .Data()[0]);
  }

  void SetLineMixing2AER(ConstVectorView d) {
    const LineShape::ModelParameters Y = {LineShape::TemperatureModel::LM_AER, d[4], d[5], d[6], d[7]};
    const LineShape::ModelParameters G = {LineShape::TemperatureModel::LM_AER, d[8], d[9], d[10], d[11]};
    for (auto& sm : mlineshapemodel.Data())
      for (Index i = 0; i < d.nelem(); i++) {
        sm.Y() = Y;
        sm.G() = G;
      }
  }

  LineShape::Model& GetLineShapeModel() { return mlineshapemodel; }
  const LineShape::Model& GetLineShapeModel() const { return mlineshapemodel; }

  void SetLineShapeModel(const LineShape::Model& lsm) { mlineshapemodel = lsm; }

  void LineShapeModelOnlyAir() {
    if (mlineshapemodel.Bath())
      mlineshapemodel = LineShape::Model(mlineshapemodel.Data().back(),
                                         mlineshapemodel.ModelType());
    else
      throw std::runtime_error(
          "Cannot keep only pure air since no air broadening exists");
  }

  void LineShapeModelRemoveSelf() {
    if (mlineshapemodel.Self()) mlineshapemodel.Remove(0);
  }

  void SetLineShapeModelParameters(const LineShape::ModelParameters& x,
                                   const String& spec,
                                   const String& var) {
    mlineshapemodel.Set(x, spec, LineShape::string2variable(var));
  }

  void SetLineShapeModelParameter(const Numeric& X,
                                  const String& spec,
                                  const String& var,
                                  const String& coeff) {
    auto x = GetLineShapeModelParameters(spec, var);
    LineShape::SingleModelParameter(x, coeff) = X;
    SetLineShapeModelParameters(x, spec, var);
  }

  LineShape::ModelParameters GetLineShapeModelParameters(
      const String& spec, const String& var) const {
    return mlineshapemodel.Get(spec, LineShape::string2variable(var));
  }

  void SetSpecial() { mstandard = false; }
  void SetStandard() { mstandard = true; }

  bool LineMixingByBand() const { return not mstandard; }

  bool LineShapeModelHasAir() const noexcept { return mlineshapemodel.Bath(); }

  void SetAirPressureBroadening(Numeric& G0,
                                Numeric& D0,
                                const Numeric& T,
                                const Numeric& P,
                                const Numeric& self_vmr) const {
    G0 = P * (mlineshapemodel.Data().front().compute(
                  T, mti0, LineShape::Variable::G0) *
                  self_vmr +
              mlineshapemodel.Data().back().compute(
                  T, mti0, LineShape::Variable::G0) *
                  (1 - self_vmr));
    D0 = P * (mlineshapemodel.Data().front().compute(
                  T, mti0, LineShape::Variable::D0) *
                  self_vmr +
              mlineshapemodel.Data().back().compute(
                  T, mti0, LineShape::Variable::D0) *
                  (1 - self_vmr));
  }

  Numeric PressureBroadeningAirBroadeningNair() const { return Nair(); }
  Numeric PressureBroadeningAirBroadeningPsf() const { return Psf(); }
  Numeric PressureBroadeningAirBroadeningAgam() const { return Agam(); }

  Numeric GetPrepInternalDerivative(const Numeric& T,
                                    const Numeric& P,
                                    const Vector& vmrs,
                                    const RetrievalQuantity& derivative) const
      noexcept {
    const bool self = derivative.Mode() == LineShape::self_broadening;
    const bool bath = derivative.Mode() == LineShape::bath_broadening;
    if (mlineshapemodel.Self() and self)
      return mlineshapemodel.GetInternalDeriv(
          T, mti0, P, 0, vmrs, derivative.PropMatType());
    else if (self) {
      const auto pos = mlineshapemodel.this_species(mqid);
      return mlineshapemodel.GetInternalDeriv(
          T, mti0, P, pos, vmrs, derivative.PropMatType());
    } else if (mlineshapemodel.Bath() and bath)
      return mlineshapemodel.GetInternalDeriv(T,
                                              mti0,
                                              P,
                                              mlineshapemodel.nelem() - 1,
                                              vmrs,
                                              derivative.PropMatType());
    else if (bath)
      return 0;
    else {
      const auto pos =
          mlineshapemodel.this_species(derivative.QuantumIdentity());
      return mlineshapemodel.GetInternalDeriv(
          T, mti0, P, pos, vmrs, derivative.PropMatType());
    }
  }

  Numeric GetInternalDerivative(const Numeric& T,
                                const Numeric& P,
                                const ConstVectorView vmrs,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const RetrievalQuantity& derivative) const {
    auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
    return GetPrepInternalDerivative(T, P, v, derivative);
  }

  LineShape::Output GetPrepShapeParams(const Numeric& T,
                                       const Numeric& P,
                                       const Vector& vmrs) const noexcept {
    auto x = mlineshapemodel.GetParams(T, mti0, P, vmrs);
    if (not do_linemixing(P)) x.Y = x.G = x.DV = 0;
    return x;
  }

  /*! Method to compute the line mixing and pressure broadening parameters */
  LineShape::Output GetShapeParams(
      const Numeric& T,
      const Numeric& P,
      const ConstVectorView vmrs,
      const ArrayOfArrayOfSpeciesTag& abs_species) const {
    auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
    return GetPrepShapeParams(T, P, v);
  }

  LineShape::Output GetPrepShapeParams_dT(const Numeric& T,
                                          const Numeric& P,
                                          const Vector& vmrs) const noexcept {
    auto x = mlineshapemodel.GetTemperatureDerivs(T, mti0, P, vmrs);
    if (not do_linemixing(P)) x.Y = x.G = x.DV = 0;
    return x;
  }

  /*! Method to compute the temperature derivatives of line mixing and pressure broadening parameters */
  LineShape::Output GetShapeParams_dT(
      const Numeric& T,
      const Numeric& P,
      const ConstVectorView vmrs,
      const ArrayOfArrayOfSpeciesTag& abs_species) const {
    auto v = mlineshapemodel.vmrs(vmrs, abs_species, mqid);
    return GetPrepShapeParams_dT(T, P, v);
  }

  /*! Method to compute the temperature derivatives of line mixing and pressure broadening parameters */
  LineShape::Output GetShapeParams_dVMR(const Numeric& T,
                                        const Numeric& P,
                                        const QuantumIdentifier& vmr_qi) const
      noexcept {
    const bool self = vmr_qi.Species() == mqid.Species();
    if (mlineshapemodel.Self() and self) {
      auto x = mlineshapemodel.GetVMRDerivs(T, mti0, P, 0);

      if (mlineshapemodel.Bath())
        x = LineShape::differenceOutput(
            x,
            mlineshapemodel.GetVMRDerivs(
                T, mti0, P, mlineshapemodel.nelem() - 1));

      if (not do_linemixing(P)) x.Y = x.G = x.DV = 0;

      return x;
    } else if (mlineshapemodel.Bath() and self)
      return {0, 0, 0, 0, 0, 0, 0, 0, 0};
    else {
      const auto pos = mlineshapemodel.this_species(vmr_qi);
      auto x = mlineshapemodel.GetVMRDerivs(T, mti0, P, pos);

      if (mlineshapemodel.Bath())
        x = LineShape::differenceOutput(
            x,
            mlineshapemodel.GetVMRDerivs(
                T, mti0, P, mlineshapemodel.nelem() - 1));

      if (not do_linemixing(P)) x.Y = x.G = x.DV = 0;
      return x;
    }
  }

  /** Zeeman Effect Data */
  void ZeemanModelInit() { mzeemanmodel = Zeeman::Model(mqid); }
  void ZeemanModelInitZero() { mzeemanmodel = Zeeman::Model(); }
  Numeric ZeemanModelStrength(Zeeman::Polarization type, Index i) const {
    return mzeemanmodel.Strength(UpperQuantumNumber(QuantumNumberType::J),
                                 LowerQuantumNumber(QuantumNumberType::J),
                                 type,
                                 i);
  }
  Numeric ZeemanModelSplitting(Zeeman::Polarization type, Index i) const {
    return mzeemanmodel.Splitting(UpperQuantumNumber(QuantumNumberType::J),
                                  LowerQuantumNumber(QuantumNumberType::J),
                                  type,
                                  i);
  }
  Index ZeemanModelLineCount(Zeeman::Polarization type) const {
    return Zeeman::nelem(UpperQuantumNumber(QuantumNumberType::J),
                         LowerQuantumNumber(QuantumNumberType::J),
                         type);
  }
  Zeeman::Model ZeemanModel() const { return mzeemanmodel; }
  Zeeman::Model& ZeemanModel() { return mzeemanmodel; }

  /** Cutoff frequency */
  const Numeric& CutOff() const { return mcutoff; }
  void SetCutOff(const Numeric& cutoff) { mcutoff = cutoff; }

  /** Line mixing pressure limit */
  const Numeric& LineMixingLimit() const { return mlinemixing_limit; }
  void SetLineMixingLimit(const Numeric& limit) { mlinemixing_limit = limit; }

  /** Line shape mirroring factor */
  const MirroringType& GetMirroringType() const { return mmirroring; }
  void SetMirroringType(const MirroringType in) { mmirroring = in; }
  String GetMirroringTypeString() const;

  /** Line shape normalization factor */
  const LineNormalizationType& GetLineNormalizationType() const {
    return mlinenorm;
  }
  void SetLineNormalizationType(const LineNormalizationType in) {
    mlinenorm = in;
  }
  String GetLineNormalizationTypeString() const;

  LineShape::Type GetLineShapeType() const noexcept {
    return mlineshapemodel.ModelType();
  }

  /** Line population type */
  const LinePopulationType& GetLinePopulationType() const {
    return mpopulation;
  }
  void SetLinePopulationType(const LinePopulationType in) { mpopulation = in; }
  String GetLinePopulationTypeString() const;

  bool IsNotSameSpecIso(const LineRecord& other) const {
    return mqid.Species() not_eq other.mqid.Species() or
           mqid.Isotopologue() not_eq other.mqid.Isotopologue();
  }

  /** Read one line from a stream associated with a HITRAN 1986-2001 file. The
    HITRAN format is as follows (directly from the HITRAN documentation):

    \verbatim
    Each line consists of 100
    bytes of ASCII text data, followed by a line feed (ASCII 10) and
    carriage return (ASCII 13) character, for a total of 102 bytes per line.
    Each line can be read using the following READ and FORMAT statement pair
    (for a FORTRAN sequential access read):

          READ(3,800) MO,ISO,V,S,R,AGAM,SGAM,E,N,d,V1,V2,Q1,Q2,IERF,IERS,
         *  IERH,IREFF,IREFS,IREFH
    800   FORMAT(I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2I3,2A9,3I1,3I2)

    Each item is defined below, with its format shown in parenthesis.

      MO  (I2)  = molecule number
      ISO (I1)  = isotopologue number (1 = most abundant, 2 = second, etc)
      V (F12.6) = frequency of transition in wavenumbers (cm-1)
      S (E10.3) = intensity in cm-1/(molec * cm-2) at 296 Kelvin
      R (E10.3) = transition probability squared in Debyes**2
      AGAM (F5.4) = air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      SGAM (F5.4) = self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      E (F10.4) = lower state energy in wavenumbers (cm-1)
      N (F4.2) = coefficient of temperature dependence of air-broadened halfwidth
      d (F8.6) = shift of transition due to pressure (cm-1)
      V1 (I3) = upper state global quanta index
      V2 (I3) = lower state global quanta index
      Q1 (A9) = upper state local quanta
      Q2 (A9) = lower state local quanta
      IERF (I1) = accuracy index for frequency reference
      IERS (I1) = accuracy index for intensity reference
      IERH (I1) = accuracy index for halfwidth reference
      IREFF (I2) = lookup index for frequency
      IREFS (I2) = lookup index for intensity
      IREFH (I2) = lookup index for halfwidth

    The molecule numbers are encoded as shown in the table below:

      0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
      6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
     12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
     18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
     24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
     30=  SF6   31=  H2S   32=HCOOH
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped. It is checked if the data record has the right
    number of characters. If not, a runtime_error is thrown.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \author Stefan Buehler */
  bool ReadFromHitran2001Stream(istream& is, const Verbosity& verbosity);

  /** LBLRTM uses the same format as HITRAN pre-2004 but also carry 
    line mixing data, so we must read it separately.

   \author Richard Larsson
   \date   2014-09-24 */
  bool ReadFromLBLRTMStream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with a HITRAN 2004 file. The
    HITRAN format is as follows:

    \verbatim
    Each line consists of 160 ASCII characters, followed by a line feed (ASCII 10)
    and carriage return (ASCII 13) character, for a total of 162 bytes per line.

    Each item is defined below, with its Fortran format shown in parenthesis.

    (I2)     molecule number
    (I1)     isotopologue number (1 = most abundant, 2 = second, etc)
    (F12.6)  vacuum wavenumbers (cm-1)
    (E10.3)  intensity in cm-1/(molec * cm-2) at 296 Kelvin
    (E10.3)  Einstein-A coefficient (s-1)
    (F5.4)   air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    (F5.4)   self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    (F10.4)  lower state energy (cm-1)
    (F4.2)   coefficient of temperature dependence of air-broadened halfwidth
    (F8.6)   air-broadened pressure shift of line transition at 296 K (cm-1)
    (A15)    upper state global quanta
    (A15)    lower state global quanta
    (A15)    upper state local quanta
    (A15)    lower state local quanta
    (I1)     uncertainty index for wavenumber
    (I1)     uncertainty index for intensity
    (I1)     uncertainty index for air-broadened half-width
    (I1)     uncertainty index for self-broadened half-width
    (I1)     uncertainty index for temperature dependence
    (I1)     uncertainty index for pressure shift
    (I2)     index for table of references correspond. to wavenumber
    (I2)     index for table of references correspond. to intensity
    (I2)     index for table of references correspond. to air-broadened half-width
    (I2)     index for table of references correspond. to self-broadened half-width
    (I2)     index for table of references correspond. to temperature dependence
    (I2)     index for table of references correspond. to pressure shift
    (A1)     flag (*) for lines supplied with line-coupling algorithm
    (F7.1)   upper state statistical weight
    (F7.1)   lower state statistical weight

     The molecule numbers are encoded as shown in the table below:

      0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=    CO
      6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=   NH3
     12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=    HI
     18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=   HCN
     24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29=  COF2
     30=  SF6   31=  H2S   32=HCOOH   33=  HO2   34=    O   35=ClONO2
     36=  NO+   37= HOBr   38= C2H4
    \endverbatim

    CH3OH is not included in ARTS because its total internal partition
    sum is not known yet.

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped. It is checked if the data record has the right number
    of characters (comment lines are ignored). If not, a runtime_error is thrown.
    If the molecule is unknown to ARTS, a warning is prompted but the
    program continues (ignoring this record). For CH3OH this
    warning will be issued even when using the regular Hitran 2004 data base
    (see above).
    If the line center is below fmin, mf is set to -1 and the caller should ignore
    this line.

    \param is Stream from which to read
    \param verbosity Verbosity
    \param fmin Skip line if mf < fmin
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \author Stefan Buehler, Hermann Berg */
  bool ReadFromHitran2004Stream(istream& is,
                                const Verbosity& verbosity,
                                const Numeric fmin = 0);

  /** Read one line from a stream associated with a MYTRAN2 file. The MYTRAN2
    format is as follows (directly taken from the abs_my.c documentation):

    \verbatim
    The MYTRAN format is as follows (FORTRAN notation):
    FORMAT(I2,I1,F13.4,1PE10.3,0P2F5.2,F10.4,2F4.2,F8.6,F6.4,2I3,2A9,4I1,3I2)
   
    Each item is defined below, with its FORMAT String shown in
    parenthesis.
   
       MO  (I2)      = molecule number
       ISO (I1)      = isotopologue number (1 = most abundant, 2 = second, etc)
    *  F (F13.4)     = frequency of transition in MHz
    *  errf (F8.4)   = error in f in MHz
       S (E10.3)     = intensity in cm-1/(molec * cm-2) at 296 K
    *  AGAM (F5.4)   = air-broadened halfwidth (HWHM) in MHz/Torr at Tref
    *  SGAM (F5.4)   = self-broadened halfwidth (HWHM) in MHz/Torr at Tref
       E (F10.4)     = lower state energy in wavenumbers (cm-1)
       N (F4.2)      = coefficient of temperature dependence of 
                       air-broadened halfwidth
    *  N_self (F4.2) = coefficient of temperature dependence of 
                       self-broadened halfwidth
    *  Tref (F7.2)   = reference temperature for AGAM and SGAM 
    *  d (F9.7)      = shift of transition due to pressure (MHz/Torr)
       V1 (I3)       = upper state global quanta index
       V2 (I3)       = lower state global quanta index
       Q1 (A9)       = upper state local quanta
       Q2 (A9)       = lower state local quanta
       IERS (I1)     = accuracy index for S
       IERH (I1)     = accuracy index for AGAM
    *  IERN (I1)     = accuracy index for N

   
    The asterisks mark entries that are different from HITRAN.

    Note that AGAM and SGAM are for the temperature Tref, while S is
    still for 296 K!
   
    The molecule numbers are encoded as shown in the table below:
   
     0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
     6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
    30=  SF6   31=  H2S   32=HCOOH   33= HO2    34=    O   35= CLONO2
    36=  NO+   37= Null   38= Null   39= Null   40=H2O_L   41= Null
    42= Null   43= OCLO   44= Null   45= Null   46=BRO     47= Null
    48= H2SO4  49=CL2O2

    All molecule numbers are from HITRAN, except for species with id's
    greater or equals 40, which are not included in HITRAN.
    (E.g.: For BrO, iso=1 is Br-79-O,iso=2 is  Br-81-O.)
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 31.10.00
    \author Axel von Engeln 
  */
  bool ReadFromMytran2Stream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with a JPL file. The JPL
    format is as follows (directly taken from the jpl documentation):

    \verbatim 
    The catalog line files are composed of 80-character lines, with one
    line entry per spectral line.  The format of each line is:

    \label{lfmt}
    \begin{tabular}{@{}lccccccccr@{}}
    FREQ, & ERR, & LGINT, & DR, & ELO, & GUP, & TAG, & QNFMT, & QN${'}$, & QN${''}$\\ 
    (F13.4, & F8.4, & F8.4, & I2, & F10.4, & I3, & I7, & I4, & 6I2, & 6I2)\\
    \end{tabular}

    \begin{tabular}{lp{4.5in}} 
    FREQ: & Frequency of the line in MHz.\\ 
    ERR: & Estimated or experimental error of FREQ in MHz.\\ 
    LGINT: &Base 10 logarithm of the integrated intensity 
    in units of \linebreak nm$^2$$\cdot$MHz at 300 K. (See Section 3 for 
    conversions to other units.)\\ 
    DR: & Degrees of freedom in the rotational partition 
    function (0 for atoms, 2 for linear molecules, and 3 for nonlinear 
    molecules).\\ 
    ELO: &Lower state energy in cm$^{-1}$ relative to the lowest energy 
    spin--rotation level in ground vibronic state.\\ 
    GUP: & Upper state degeneracy.\\ 
    TAG: & Species tag or molecular identifier. 
    A negative value flags that the line frequency has 
    been measured in the laboratory.  The absolute value of TAG is then the 
    species tag and ERR is the reported experimental error.  The three most 
    significant digits of the species tag are coded as the mass number of the 
    species, as explained above.\\ 
    QNFMT: &Identifies the format of the quantum numbers 
    given in the field QN. These quantum number formats are given in Section 5 
    and are different from those in the first two editions of the catalog.\\ 
    QN${'}$: & Quantum numbers for the upper state coded 
    according to QNFMT.\\ 
    QN${''}$: & Quantum numbers for the lower state.\\
    \end{tabular} 
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped (unused in jpl).

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 01.11.00
    \author Axel von Engeln */
  bool ReadFromJplStream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with an ARTSCAT-3 file.

      Format: see Documentation of class LineRecord

      The function attempts to read a line of data from the
      catalogue. It returns false if it succeeds. Otherwise, if eof is
      reached, it returns true. If an error occurs, a runtime_error is
      thrown. When the function looks for a data line, comment lines are
      automatically skipped.

      \param is Stream from which to read
      \exception runtime_error Some error occured during the read
      \return false=ok (data returned), true=eof (no data returned)

      \date   2000-12-15
      \author Oliver Lemke

      
      \date   2001-06-20
      \author Stefan Buehler
*/
  bool ReadFromArtscat3Stream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with an ARTSCAT-4 file.
   
   Format: see Documentation of class LineRecord
   
   The function attempts to read a line of data from the
   catalogue. It returns false if it succeeds. Otherwise, if eof is
   reached, it returns true. If an error occurs, a runtime_error is
   thrown. When the function looks for a data line, comment lines are
   automatically skipped.
   
   \param is Stream from which to read
   \exception runtime_error Some error occured during the read
   \return false=ok (data returned), true=eof (no data returned)
   
   \date   2012-02-10
   \author Oliver Lemke
   */
  bool ReadFromArtscat4Stream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with an ARTSCAT-5 file.
   
   Format: see Documentation of class LineRecord
   
   The function attempts to read a line of data from the
   catalogue. It returns false if it succeeds. Otherwise, if eof is
   reached, it returns true. If an error occurs, a runtime_error is
   thrown. When the function looks for a data line, comment lines are
   automatically skipped.
   
   \param is Stream from which to read
   \exception runtime_error Some error occured during the read
   \return false=ok (data returned), true=eof (no data returned)
   
   \date   2014-12-01
   \author Oliver Lemke
   */
  bool ReadFromArtscat5Stream(istream& is, const Verbosity& verbosity);

 private:
  // Version number:
  Index mversion;

  // Identity of line
  QuantumIdentifier mqid;

  // The line center frequency in Hz:
  Numeric mf;

  // The line intensity in m^2/Hz:
  Numeric mi0;

  // Reference temperature for I0 in K:
  Numeric mti0;

  // Lower state energy in Joules:
  Numeric melow;

  // Lower state vibrational energy in Joules:
  Numeric mevlow;  // NOTE: Not stored in binary data

  // Lower state vibrational energy index:
  Index mnlte_lower_index;  // NOTE: Not stored in binary data

  // Upper state vibrational energy in Joules:
  Numeric mevupp;  // NOTE: Not stored in binary data

  // Upper state vibrational energy index:
  Index mnlte_upper_index;  // NOTE: Not stored in binary data

  // Einstein A-coefficient in 1/s:
  Numeric ma;

  // Upper state stat. weight:
  Numeric mgupper;

  // Lower state stat. weight:
  Numeric mglower;

  /** Line function data (pressure broadening and line mixing)
   *
   *  FIXME Richard, please explain.
   * */
  LineShape::Model mlineshapemodel;

  /** FIXME Richard, please explain.   */
  bool mstandard;

  /** Zeeman effect model class */
  Zeeman::Model mzeemanmodel;

  /** Cutoff frequency */
  Numeric mcutoff;

  /** Line mixing pressure limit */
  Numeric mlinemixing_limit;

  /** Line shape mirroring effect type */
  MirroringType mmirroring;

  /** Line shape normalization type */
  LineNormalizationType mlinenorm;

  /** Line LTE/NLTE type */
  LinePopulationType mpopulation;
};

/** Output operator for LineRecord. The result should look like a
    catalogue line.

    \author Stefan Buehler */
ostream& operator<<(ostream& os, const LineRecord& lr);

//======================================================================
//         Typedefs for LineRecord Arrays
//======================================================================

/** Holds a list of spectral line data.
    \author Stefan Buehler */
typedef Array<LineRecord> ArrayOfLineRecord;

/** Holds a lists of spectral line data for each tag group.
    Dimensions: (tag_groups.nelem()) (number of lines for this tag)
    \author Stefan Buehler */
typedef Array<Array<LineRecord> > ArrayOfArrayOfLineRecord;

LinePopulationType LinePopulationTypeFromString(const String& in);

MirroringType MirroringTypeFromString(const String& in);

LineNormalizationType LineNormalizationTypeFromString(const String& in);

#endif  // linerecord_h
