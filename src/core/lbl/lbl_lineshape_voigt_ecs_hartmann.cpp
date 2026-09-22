#include "lbl_lineshape_voigt_ecs_hartmann.h"

#include <arts_conversions.h>
#include <atm.h>
#include <wigner_functions.h>

#include <cmath>

namespace lbl::voigt::ecs::hartmann {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

namespace {
void validate_rotational_state(const Rational J, const Rational l) {
  ARTS_USER_ERROR_IF(J.denom != 1 or l.denom != 1 or J < 0 or abs(l) > J,
                     "Hartmann ECS requires integer J >= 0 and integer |l| <= J; got J={}, l={}",
                     J,
                     l);
}

Numeric wig3(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER3(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER6(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

std::function<Numeric(Rational)> erot_selection(const SpeciesIsotope& isot) {
  if (isot == "CO2-626"_isot) {
    return [](const Rational J) -> Numeric { return Conversion::kaycm2joule(0.39021) * Numeric(J * (J + 1)); };
  }

  ARTS_USER_ERROR(
      "Hartmann ECS currently supports rotational energies only for CO2-626; {} requires "
      "isotopologue-specific rotational energies and validated collision parameters.",
      isot.FullName())
  return [](const Rational J) -> Numeric { return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN(); };
}
}  // namespace

Numeric reduced_dipole(const Rational Jf, const Rational Ji, const Rational lf, const Rational li, const Rational k) {
  validate_rotational_state(Jf, lf);
  validate_rotational_state(Ji, li);
  if (not iseven(Jf + lf + 1)) return -sqrtr(2 * Jf + 1) * wigner3j(Jf, k, Ji, lf, li - lf, -li);
  return +sqrtr(2 * Jf + 1) * wigner3j(Jf, k, Ji, lf, li - lf, -li);
}

void relaxation_matrix_offdiagonal(MatrixView&                     W,
                                   const QuantumIdentifier&        bnd_qid,
                                   const band_data&                bnd,
                                   const ArrayOfIndex&             sorting,
                                   const SpeciesEnum               broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector&                   dipr,
                                   const AtmPoint&                 atm) {
  using Conversion::kelvin2joule;

  const Size n = bnd.size();
  if (not n) return;

  // This linear-rotor kernel has no elastic L=0 basis rate.  Distinct
  // transitions with the same rotational pair require additional state labels
  // and collision dynamics; evaluating the power law at L=0 is singular.
  for (Size i = 0; i < n; ++i) {
    const auto& Ji = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J);
    for (Size j = i + 1; j < n; ++j) {
      const auto& Jj = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J);
      ARTS_USER_ERROR_IF(Ji.upper == Jj.upper and Ji.lower == Jj.lower,
                         "Hartmann ECS does not support distinct lines with the same rotational pair "
                         "(J upper={}, J lower={}); an elastic L=0 collision model is required.",
                         Ji.upper,
                         Ji.lower)
    }
  }

  // These are constant for a band
  auto&    l2 = bnd_qid.state.at(QuantumNumberType::l2);
  Rational li = l2.upper;
  Rational lf = l2.lower;

  for (const auto& ln : bnd) {
    const auto& J = ln.qn.at(QuantumNumberType::J);
    validate_rotational_state(J.upper, li);
    validate_rotational_state(J.lower, lf);
  }

  using std::swap;
  const bool swap_order = li > lf;
  if (swap_order) swap(li, lf);
  if (abs(li - lf) > 1) return;

  const Numeric T = atm.temperature;

  const auto erot = erot_selection(bnd_qid.isot);

  const std::array rats{bnd.max(QuantumNumberType::J), li, lf};
  const int        maxL = wigner_init_size(rats);

  const auto Om = [&]() {
    Vector out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Omega(
          atm.temperature,
          bnd.front().ls.T0,
          broadening_species == SpeciesEnum::Bath ? atm.mean_mass() : atm.mean_mass(broadening_species),
          bnd_qid.isot.mass,
          erot(Rational{i}),
          erot(Rational{i - 2}));
    return out;
  }();

  const auto Q = [&]() {
    Vector out(maxL, 0.0);
    // Only positive transfer channels are defined by the power law.
    for (Index i = 1; i < maxL; i++)
      out[i] = rovib_data.Q(Rational{i}, atm.temperature, bnd.front().ls.T0, erot(Rational{i}));
    return out;
  }();

  for (Index L = 0; L < maxL; ++L) {
    ARTS_USER_ERROR_IF(not std::isfinite(Om[L]) or Om[L] <= 0 or not std::isfinite(Q[L]),
                       "Invalid ECS basis rate or adiabaticity factor at L={} for {} and {}",
                       L,
                       bnd_qid.isot,
                       broadening_species)
  }

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    auto&    J  = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J);
    Rational Ji = J.upper;
    Rational Jf = J.lower;
    if (swap_order) swap(Ji, Jf);

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;
      auto&    J_p  = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J);
      Rational Ji_p = J_p.upper;
      Rational Jf_p = J_p.lower;
      if (swap_order) swap(Ji_p, Jf_p);

      // Select upper quantum number
      if (Jf_p > Jf) continue;

      Index L         = std::max(std::abs((Ji - Ji_p).toIndex()), std::abs((Jf - Jf_p).toIndex()));
      L              += L % 2;
      const Index Lf  = std::min((Ji + Ji_p).toIndex(), (Jf + Jf_p).toIndex());

      Numeric sum = 0;
      for (; L <= Lf; L += 2) {
        const Numeric a  = wig3(Ji, Ji_p, Rational{L}, li, -li, Rational{0});
        const Numeric b  = wig3(Jf, Jf_p, Rational{L}, lf, -lf, Rational{0});
        const Numeric c  = wig6(Ji, Jf, Rational{1}, Jf_p, Ji_p, Rational{L});
        sum             += a * b * c * Numeric(2 * L + 1) * Q[L] / Om[L];
      }
      const Numeric ECS  = Om[Ji.toIndex()];
      const Numeric scl  = ECS * Numeric(2 * Ji_p + 1) * sqrtr((2 * Jf + 1) * (2 * Jf_p + 1));
      sum               *= scl;

      // Add to W and rescale to upwards element by the populations
      W[j, i] = sum;
      W[i, j] = sum * std::exp((erot(Jf_p) - erot(Jf)) / kelvin2joule(T));
    }
  }
  arts_wigner_thread_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")

  // The sequential correction retains the historical truncated-band closure.
  // In particular it cannot enforce the final column's sum rule.  Do not hide
  // overflow or invalid rates by treating a non-finite denominator as zero.
  for (Size i = 0; i < n; ++i) {
    ARTS_USER_ERROR_IF(not std::isfinite(dipr[i]), "Non-finite ECS reduced dipole for line {}", sorting[i])
    for (Size j = 0; j < n; ++j) {
      ARTS_USER_ERROR_IF(not std::isfinite(W[j, i]),
                         "Non-finite ECS relaxation matrix element ({}, {}) before sum-rule correction",
                         j,
                         i)
    }
  }

  // Sum rule correction
  for (Size i = 0; i < n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;

    for (Size j = 0; j < n; j++) {
      if (j > i) {
        sumlw += dipr[j] * W[j, i];
      } else {
        sumup += dipr[j] * W[j, i];
      }
    }

    const Rational Ji = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J).lower;
    ARTS_USER_ERROR_IF(
        not std::isfinite(sumlw) or not std::isfinite(sumup), "Non-finite ECS sum-rule sums for line {}", sorting[i])
    ARTS_USER_ERROR_IF(sumlw != 0 and not std::isfinite(-sumup / sumlw),
                       "ECS sum-rule correction overflows for line {}; the supplied band and widths "
                       "do not define a stable correction.",
                       sorting[i])

    for (Size j = i + 1; j < n; j++) {
      const Rational Jj = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J).lower;
      if (sumlw == 0) {
        W[j, i] = 0.0;
        W[i, j] = 0.0;
      } else {
        W[j, i] *= -sumup / sumlw;
        W[i, j]  = W[j, i] * std::exp((erot(Ji) - erot(Jj)) / kelvin2joule(T));  // This gives LTE
      }
    }
  }

  for (Size i = 0; i < n; ++i) {
    for (Size j = 0; j < n; ++j) {
      ARTS_USER_ERROR_IF(not std::isfinite(W[j, i]),
                         "Non-finite ECS relaxation matrix element ({}, {}) after sum-rule correction",
                         j,
                         i)
    }
  }
}
}  // namespace lbl::voigt::ecs::hartmann
