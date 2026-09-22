#include "lbl_lineshape_voigt_ecs_makarov.h"

#include <arts_conversions.h>
#include <wigner_functions.h>

#include <cmath>

namespace lbl::voigt::ecs::makarov {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

namespace {
void validate_rotational_pair(const Rational Ju, const Rational Jl, const Rational N) {
  ARTS_USER_ERROR_IF(
      N.denom != 1 or N <= 0 or iseven(N), "Makarov O2-66 ECS requires a positive odd integer N, got {}", N);
  ARTS_USER_ERROR_IF(Ju.denom != 1 or Jl.denom != 1 or Ju < 0 or Jl < 0 or abs(Ju - N) > 1 or abs(Jl - N) > 1,
                     "Makarov O2-66 ECS requires integer J >= 0 with |J-N| <= 1; got Ju={}, Jl={}, N={}",
                     Ju,
                     Jl,
                     N);
  ARTS_USER_ERROR_IF(abs(Ju - Jl) > 1 or (Ju == 0 and Jl == 0),
                     "Makarov O2-66 ECS requires a rank-1 transition with |Ju-Jl| <= 1 and excludes J=0 to J=0; "
                     "got Ju={}, Jl={}",
                     Ju,
                     Jl);
}

Numeric wig3(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER3(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER6(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}
}  // namespace

Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N) {
  validate_rotational_pair(Ju, Jl, N);
  return (iseven(Jl + N) ? 1 : -1) * sqrtr(6 * (2 * Jl + 1) * (2 * Ju + 1)) *
         wigner6j(Rational{1}, Rational{1}, Rational{1}, Jl, Ju, N);
};

namespace {
/*! Compute the rotational energy of ground-state O2 at N and J
 * 
 * If the template argument evaluates true, the erot<false>(1, 0)
 * energy is removed from the output of erot<false>(N, J).
 * 
 * @param[in] N Main rotational number
 * @param[in] j Main rotational number plus spin (if j < 0 then J=N)
 * @return Rotational energy in Joule
 */
template <bool rescale_pure_rotational = true>
constexpr Numeric erot(const Rational N, const Rational j = Rational{-1}) try {
  const Rational J = j < 0 ? N : j;

  if constexpr (rescale_pure_rotational) {
    return erot<false>(N, J) - erot<false>(Rational{1}, Rational{0});
  } else {
    using Conversion::mhz2joule;
    using Math::pow2;
    using Math::pow3;

    constexpr Numeric B0  = 43100.4425e0;
    constexpr Numeric D0  = .145123e0;
    constexpr Numeric H0  = 3.8e-08;
    constexpr Numeric xl0 = 59501.3435e0;
    constexpr Numeric xg0 = -252.58633e0;
    constexpr Numeric xl1 = 0.058369e0;
    constexpr Numeric xl2 = 2.899e-07;
    constexpr Numeric xg1 = -2.4344e-04;
    constexpr Numeric xg2 = -1.45e-09;

    const Numeric XN      = static_cast<Numeric>(N);
    const Numeric XX      = XN * (XN + 1);
    const Numeric xlambda = xl0 + xl1 * XX + xl2 * pow2(XX);
    const Numeric xgama   = xg0 + xg1 * XX + xg2 * pow2(XX);
    const Numeric C1      = B0 * XX - D0 * pow2(XX) + H0 * pow3(XX);

    if (J < N) {
      if (N == 1)  // erot<false>(1, 0)
        return mhz2joule(C1 - (xlambda + B0 * (2. * XN - 1.) + xgama * XN));
      return mhz2joule(C1 - (xlambda + B0 * (2. * XN - 1.) + xgama * XN) +
                       std::sqrt(pow2(B0 * (2. * XN - 1.)) + pow2(xlambda) - 2. * B0 * xlambda));
    }
    if (J > N)
      return mhz2joule(C1 - (xlambda - B0 * (2. * XN + 3.) - xgama * (XN + 1.)) -
                       std::sqrt(pow2(B0 * (2. * XN + 3.)) + pow2(xlambda) - 2. * B0 * xlambda));
    return mhz2joule(C1);
  }
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd) {
  ARTS_USER_ERROR_IF(bnd_qid.isot != "O2-66"_isot, "Makarov ECS currently supports only O2-66, got {}", bnd_qid.isot);
  const auto& S = bnd_qid.state.at(QuantumNumberType::S);
  ARTS_USER_ERROR_IF(S.upper != 1 or S.lower != 1,
                     "Makarov O2-66 ECS requires electron spin S=1 in both states; got upper={}, lower={}",
                     S.upper,
                     S.lower);
  for (const auto& ln : bnd) {
    const auto& J = ln.qn.at(QuantumNumberType::J);
    const auto& N = ln.qn.at(QuantumNumberType::N);
    ARTS_USER_ERROR_IF(N.upper != N.lower,
                       "Makarov O2-66 ECS implements microwave transitions with unchanged N; got upper={}, lower={}",
                       N.upper,
                       N.lower);
    validate_rotational_pair(J.upper, J.lower, N.upper);
  }
}

void relaxation_matrix_offdiagonal(MatrixView&                     W,
                                   const QuantumIdentifier&        bnd_qid,
                                   const band_data&                bnd,
                                   const ArrayOfIndex&             sorting,
                                   const SpeciesEnum               broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector&                   dipr,
                                   const AtmPoint&                 atm) try {
  using Conversion::kelvin2joule;

  if (bnd.size() == 0) return;
  validate_band(bnd_qid, bnd);

  const auto bk = [](const Rational& r) -> Numeric { return sqrtr(2 * r + 1); };

  const auto n = bnd.size();

  auto&          S  = bnd_qid.state.at(QuantumNumberType::S);
  const Rational Si = S.upper;
  const Rational Sf = S.lower;

  const std::array rats{bnd.max(QuantumNumberType::J), bnd.max(QuantumNumberType::N), Si, Sf};
  const int        maxL = wigner_init_size(rats);

  const auto Om = [&, maxL]() {
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

  const auto Q = [&, maxL]() {
    Vector out(maxL, 0.0);
    // The Makarov sum starts at L=2; the power law is undefined at L=0.
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
    auto& J = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J);
    auto& N = bnd.lines[sorting[i]].qn.at(QuantumNumberType::N);

    const Rational Ji = J.upper;
    const Rational Jf = J.lower;
    const Rational Ni = N.upper;
    const Rational Nf = N.lower;

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;

      auto& J_p = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J);
      auto& N_p = bnd.lines[sorting[j]].qn.at(QuantumNumberType::N);

      const Rational Ji_p = J_p.upper;
      const Rational Jf_p = J_p.lower;
      const Rational Ni_p = N_p.upper;
      const Rational Nf_p = N_p.lower;

      if (Jf_p > Jf) continue;

      // Tran etal 2006 symbol with modifications:
      //    1) [Ji] * [Ji_p] instead of [Ji_p] ^ 2 in partial accordance with Makarov etal 2013
      Numeric       sum = 0;
      const Numeric scl = (iseven(Ji_p + Ji + 1) ? 1 : -1) * bk(Ni) * bk(Nf) * bk(Nf_p) * bk(Ni_p) * bk(Jf) * bk(Jf_p) *
                          bk(Ji) * bk(Ji_p);
      const auto [L0, L1] =
          wigner_limits(wigner3j_limits<3>(Ni_p, Ni), {Rational(2), Rational{std::numeric_limits<Index>::max()}});
      for (Rational L = L0; L <= L1; L += 2) {
        const Numeric a  = wig3(Ni_p, Ni, L, Rational{0}, Rational{0}, Rational{0});
        const Numeric b  = wig3(Nf_p, Nf, L, Rational{0}, Rational{0}, Rational{0});
        const Numeric c  = wig6(L, Ji, Ji_p, Si, Ni_p, Ni);
        const Numeric d  = wig6(L, Jf, Jf_p, Sf, Nf_p, Nf);
        const Numeric e  = wig6(L, Ji, Ji_p, Rational{1}, Jf_p, Jf);
        sum             += a * b * c * d * e * Numeric(2 * L + 1) * Q[L.toIndex()] / Om[L.toIndex()];
      }
      sum *= scl * Om[Ni.toIndex()];

      // Add to W and rescale to upwards element by the populations
      W[i, j] = sum;
      W[j, i] = sum * std::exp((bnd.lines[sorting[j]].e0 - bnd.lines[sorting[i]].e0) / kelvin2joule(atm.temperature));
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

    ARTS_USER_ERROR_IF(
        not std::isfinite(sumlw) or not std::isfinite(sumup), "Non-finite ECS sum-rule sums for line {}", sorting[i])
    ARTS_USER_ERROR_IF(sumlw != 0 and not std::isfinite(-sumup / sumlw),
                       "ECS sum-rule correction overflows for line {}; the supplied band and widths "
                       "do not define a stable correction.",
                       sorting[i])

    for (Size j = i + 1; j < n; j++) {
      if (sumlw == 0) {
        W[j, i] = 0.0;
        W[i, j] = 0.0;
      } else {
        W[j, i] *= -sumup / sumlw;
        W[i, j] =
            W[j, i] * std::exp((bnd.lines[sorting[i]].e0 - bnd.lines[sorting[j]].e0) / kelvin2joule(atm.temperature));
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
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs::makarov