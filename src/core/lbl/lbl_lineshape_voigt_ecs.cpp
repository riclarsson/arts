#include "lbl_lineshape_voigt_ecs.h"

#include <arts_omp.h>
#include <atm.h>
#include <configtypes.h>
#include <debug.h>
#include <isotopologues.h>
#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <ranges>
#include <stdexcept>

#include "lbl_lineshape_linemixing.h"
#include "lbl_lineshape_model.h"
#include "lbl_lineshape_voigt_ecs_hartmann.h"
#include "lbl_lineshape_voigt_ecs_makarov.h"

#undef WIGNER3
#undef WIGNER6

namespace lbl::voigt::ecs {
ComputeData::ComputeData(const ConstVectorView&   f_grid,
                         const AtmPoint&          atm,
                         const Vector2&           los,
                         const ZeemanPolarization pol)
    : scl(f_grid.size()), shape(f_grid.size()) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 scl.begin(),
                 [N = number_density(atm.pressure, atm.temperature), T = atm.temperature](auto f) {
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -N * f * std::expm1(-r);
                 });

  update_zeeman(los, atm.mag, pol);
}

void ComputeData::update_zeeman(const Vector2& los, const Vector3& mag, const ZeemanPolarization pol) {
  npm = zeeman::norm_view(pol, mag, los);
}

void ComputeData::core_calc_eqv() {
  // The kernels use a row-rate convention; the spectral operator is Ws^T.
  // Preserve Ws and Vs so evaluation is repeatable after a single adaptation.
  // A resolvent derivative can use dR = R (dM) R without differentiating V.
  const auto n = pop.size();
  const auto m = vmrs.size();
  ARTS_USER_ERROR_IF(dip.size() != n or Ws.npages() != static_cast<Index>(m) or Ws.nrows() != static_cast<Index>(n) or
                         Ws.ncols() != static_cast<Index>(n),
                     "Inconsistent ECS matrix and population dimensions")
  eqv_strs.resize(m, n);
  eqv_vals.resize(m, n);
  Vs.resize(m, n, n);
  eigenvector_rcond.resize(m);
  eqv_strs          = 0;
  eigenvector_rcond = 1;
  if (n == 0) return;

  ComplexVector rhs(n), coefficients(n);
  for (Size j = 0; j < n; ++j) { rhs[j] = pop[j] * dip[j]; }
  complex_diagonalize_workdata workspace(n);

  for (Size k = 0; k < m; ++k) {
    auto          V = Vs[k];
    ComplexMatrix W{transpose(Ws[k])};
    auto          eqv_str = eqv_strs[k];
    auto          eqv_val = eqv_vals[k];

    // Remove the optical carrier before solving for much smaller shifts/widths.
    const Numeric center = W[0, 0].real();
    for (Size j = 0; j < n; ++j) { W[j, j] -= center; }
    Numeric matrix_norm = 0;
    for (Size i = 0; i < n; ++i) {
      Numeric row_sum = 0;
      for (Size j = 0; j < n; ++j) { row_sum += std::abs(W[i, j]); }
      matrix_norm = std::max(matrix_norm, row_sum);
    }
    ARTS_USER_ERROR_IF(not std::isfinite(matrix_norm), "ECS centered matrix norm overflowed")

    diagonalize(V, eqv_val, W, workspace);
    // Nearly defective modes have unstable individual residues even when their
    // sum is finite. Do not silently return spectra from such a decomposition.
    eigenvector_rcond[k]            = solve(coefficients, V, rhs, 1e-12);
    const Numeric damping_tolerance = 64 * std::numeric_limits<Numeric>::epsilon() * matrix_norm;
    for (Size i = 0; i < n; ++i) {
      ARTS_USER_ERROR_IF(eqv_val[i].imag() < -damping_tolerance,
                         "ECS relaxation matrix has a negative damping eigenvalue: {} Hz",
                         eqv_val[i])
      if (eqv_val[i].imag() < 0) eqv_val[i].imag(0);
      eqv_val[i]         += center;
      Complex projection  = 0;
      for (Size j = 0; j < n; ++j) { projection += dip[j] * V[j, i]; }
      eqv_str[i] = projection * coefficients[i];
      ARTS_USER_ERROR_IF(not std::isfinite(eqv_str[i].real()) or not std::isfinite(eqv_str[i].imag()),
                         "Non-finite ECS equivalent-line strength")
    }
  }
}

void ComputeData::core_calc(const ConstVectorView& f_grid) try {
  core_calc_eqv();

  const auto m = vmrs.size();
  const auto n = f_grid.size();
  shape.resize(n);
  shape = 0;

  for (Size k = 0; k < m; k++) {
    if (vmrs[k] == 0) continue;
    for (Size i = 0; i < eqv_strs[k].size(); i++) {
      const Numeric gamd = gd_fac * eqv_vals[k][i].real();
      ARTS_USER_ERROR_IF(not std::isfinite(gamd) or gamd <= 0,
                         "ECS Gaussian 1/e half-width must be positive and finite, got {} Hz",
                         gamd)
      const Numeric cte = 1 / gamd;
      for (Size iv = 0; iv < n; iv++) {
        const Complex z  = (eqv_vals[k][i] - f_grid[iv]) * cte;
        shape[iv]       += vmrs[k] * eqv_strs[k][i] * Faddeeva::w(z) / gamd;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
void get_vmrs(VectorView vmrs, const line_shape::model::map_t& mod, const AtmPoint& atm) {
  std::transform(mod.begin(), mod.end(), vmrs.begin(), [&atm](const auto& m) {
    const Numeric x = m.first == SpeciesEnum::Bath or not atm.has(m.first) ? 0.0 : atm[m.first];
    ARTS_USER_ERROR_IF(not std::isfinite(x) or x < 0, "Invalid ECS broadener VMR for {}: {}", m.first, x)
    return x;
  });

  const Numeric total = sum(vmrs);
  ARTS_USER_ERROR_IF(not std::isfinite(total), "Non-finite total ECS broadener VMR")
  const Size bath_spec = std::distance(mod.begin(), mod.find(SpeciesEnum::Bath));
  if (bath_spec != mod.size()) {
    ARTS_USER_ERROR_IF(total > 1 + 1e-12, "ECS explicit broadener VMRs exceed unity: {}", total)
    vmrs[bath_spec] = std::max(Numeric{0}, 1 - total);
  } else if (total > 0) {
    vmrs /= total;
  }
  // With no bath and no present collision partners the collision matrix is zero,
  // matching ordinary line-shape mixing, rather than producing 0/0.
}

using Offdiagonal = decltype(&hartmann::relaxation_matrix_offdiagonal);
Offdiagonal offdiagonal_kernel(LineByLineLineshape lineshape) {
  using enum LineByLineLineshape;
  switch (lineshape) {
    case VP_ECS_MAKAROV:  return makarov::relaxation_matrix_offdiagonal;
    case VP_ECS_HARTMANN: return hartmann::relaxation_matrix_offdiagonal;
    default:              ARTS_USER_ERROR("Unknown ECS line shape {}", lineshape)
  }
}

Numeric reduced_dipole(const QuantumIdentifier& qid, const band_data& bnd, const line& ln) {
  const auto& J = ln.qn.at(QuantumNumberType::J);
  using enum LineByLineLineshape;
  switch (bnd.lineshape) {
    case VP_ECS_MAKAROV:  return makarov::reduced_dipole(J.upper, J.lower, ln.qn.at(QuantumNumberType::N).upper);
    case VP_ECS_HARTMANN: {
      const auto& l = qid.state.at(QuantumNumberType::l2);
      return hartmann::reduced_dipole(J.upper, J.lower, l.upper, l.lower);
    }
    default: ARTS_USER_ERROR("Unknown ECS line shape {}", bnd.lineshape)
  }
}

Numeric closure_residual(ConstMatrixView W, ConstVectorView d) {
  Numeric result = 0;
  for (Index i = 0; i < W.ncols(); ++i) {
    Numeric residual = 0, scale = 0;
    for (Index j = 0; j < W.nrows(); ++j) {
      const Numeric term = d[j] * W[j, i];
      ARTS_USER_ERROR_IF(not std::isfinite(term), "Non-finite ECS optical sum-rule term")
      residual += term;
      scale    += std::abs(term);
    }
    ARTS_USER_ERROR_IF(not std::isfinite(residual) or not std::isfinite(scale),
                       "Non-finite ECS optical sum-rule residual")
    if (scale > 0) result = std::max(result, std::abs(residual) / scale);
  }
  return result;
}
}  // namespace

void ComputeData::adapt_multi(const QuantumIdentifier&        bnd_qid,
                              const band_data&                bnd,
                              const LinemixingSpeciesEcsData& rovib_data,
                              const AtmPoint&                 atm,
                              const bool                      presorted) {
  adapt(bnd_qid, bnd, rovib_data, atm, presorted, true);
}

void ComputeData::adapt_single(const QuantumIdentifier&        bnd_qid,
                               const band_data&                bnd,
                               const LinemixingSpeciesEcsData& rovib_data,
                               const AtmPoint&                 atm,
                               const bool                      presorted) {
  adapt(bnd_qid, bnd, rovib_data, atm, presorted, false);
}

void ComputeData::adapt(const QuantumIdentifier&        bnd_qid,
                        const band_data&                bnd,
                        const LinemixingSpeciesEcsData& rovib_data,
                        const AtmPoint&                 atm,
                        const bool                      presorted,
                        const bool                      per_broadener) try {
  const Size n = bnd.size();
  ARTS_USER_ERROR_IF(n == 0, "Cannot adapt an empty ECS band")
  ARTS_USER_ERROR_IF(
      not std::isfinite(atm.temperature) or atm.temperature <= 0 or not std::isfinite(atm.pressure) or atm.pressure < 0,
      "ECS requires positive finite temperature and nonnegative finite pressure")
  const auto& models          = bnd.front().ls.single_models;
  const Size  broadener_count = models.size();
  ARTS_USER_ERROR_IF(broadener_count == 0, "No broadening species in the ECS band")
  const Size m      = per_broadener ? broadener_count : 1;
  const auto kernel = offdiagonal_kernel(bnd.lineshape);
  ARTS_USER_ERROR_IF(bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN and bnd_qid.isot != "CO2-626"_isot,
                     "Hartmann ECS currently supports rotational energies only for CO2-626, got {}",
                     bnd_qid.isot)
  ARTS_USER_ERROR_IF(bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV and bnd_qid.isot != "O2-66"_isot,
                     "Makarov ECS currently supports only the O2-66 microwave band, got {}",
                     bnd_qid.isot)
  if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) makarov::validate_band(bnd_qid, bnd);

  for (const auto& ln : bnd) {
    ARTS_USER_ERROR_IF(
        ln.ls.single_models.size() != broadener_count or
            not stdr::all_of(models | stdv::keys, [&](auto spec) { return ln.ls.single_models.contains(spec); }),
        "All lines in an ECS band must have the same broadening species")
    ARTS_USER_ERROR_IF(not std::isfinite(ln.ls.T0) or ln.ls.T0 <= 0 or ln.ls.T0 != bnd.front().ls.T0,
                       "All lines in an ECS band must have the same positive reference temperature")
    ARTS_USER_ERROR_IF(not std::isfinite(ln.f0) or ln.f0 <= 0 or not std::isfinite(ln.a) or ln.a < 0 or
                           not std::isfinite(ln.gu) or ln.gu <= 0 or not std::isfinite(ln.e0),
                       "Invalid ECS line frequency, Einstein A, statistical weight, or lower-state energy")
  }
  if (presorted) {
    ARTS_USER_ERROR_IF(sort.size() != n, "ECS presorting requires a previous adaptation of the same band size")
    auto indices = sort;
    stdr::sort(indices);
    for (Size i = 0; i < n; ++i) {
      ARTS_USER_ERROR_IF(indices[i] != static_cast<Index>(i), "Invalid ECS sorting permutation")
    }
  }

  pop.resize(n);
  dip.resize(n);
  dipr.resize(n);
  sort.resize(n);
  Wimag.resize(n, n);
  vmrs.resize(m);
  eqv_strs.resize(m, n);
  eqv_vals.resize(m, n);
  Ws.resize(m, n, n);
  Vs.resize(m, n, n);
  sum_rule_residual.resize(broadener_count);
  Ws                = 0;
  eqv_strs          = 0;
  eqv_vals          = 0;
  sum_rule_residual = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared * atm.temperature / bnd_qid.isot.mass);
  ARTS_USER_ERROR_IF(not std::isfinite(gd_fac) or gd_fac <= 0, "Invalid ECS Doppler width factor")
  const Numeric QT = PartitionFunctions::Q(atm.temperature, bnd_qid.isot);
  ARTS_USER_ERROR_IF(not std::isfinite(QT) or QT <= 0, "Invalid ECS partition function: {}", QT)
  for (Size i = 0; i < n; ++i) {
    const auto& ln = bnd.lines[i];
    pop[i]         = ln.gu * std::exp(-ln.e0 / (Constant::k * atm.temperature)) / QT;
    dipr[i]        = reduced_dipole(bnd_qid, bnd, ln);
    ARTS_USER_ERROR_IF(not std::isfinite(dipr[i]), "Non-finite ECS reduced dipole")
    dip[i] = std::copysign(0.5 * Constant::c * std::sqrt(ln.a / (Math::pow3(ln.f0) * Constant::two_pi)), dipr[i]);
    ARTS_USER_ERROR_IF(
        not std::isfinite(pop[i]) or not std::isfinite(dip[i]) or not std::isfinite(ln.f0 * pop[i] * dip[i] * dip[i]),
        "Non-finite ECS population, dipole, or sorting strength")
  }

  if (not presorted) {
    stdr::iota(sort, 0);
    stdr::sort(stdv::zip(sort, pop, dip, dipr), stdr::greater(), [&](const auto& v) {
      const auto& [i, pop_i, dip_i, dipr_i] = v;
      return bnd.lines[i].f0 * pop_i * dip_i * dip_i;
    });
  } else {
    const auto reorder = [this](const Vector& vec) {
      Vector out(vec.size());
      for (Size i = 0; i < sort.size(); ++i) { out[i] = vec[sort[i]]; }
      return out;
    };
    pop  = reorder(pop);
    dip  = reorder(dip);
    dipr = reorder(dipr);
  }

  Vector fractions(broadener_count);
  get_vmrs(fractions, models, atm);
  if (per_broadener)
    vmrs = fractions;
  else
    vmrs = 1;

  Size i = 0;
  for (auto spec : models | stdv::keys) {
    const Size    page   = per_broadener ? i : 0;
    const Numeric weight = per_broadener ? 1 : fractions[i];
    if (weight == 0) {
      ++i;
      continue;
    }
    const auto data = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(data == rovib_data.end(), "No ECS collision data for species {}", spec)
    Wimag = 0;
    for (Size k = 0; k < n; ++k) {
      const auto&   ln    = bnd.lines[sort[k]];
      const auto&   model = ln.ls.single_models.at(spec);
      const Numeric width = model.G0(ln.ls.T0, atm.temperature, atm.pressure);
      const Numeric shift = model.D0(ln.ls.T0, atm.temperature, atm.pressure);
      ARTS_USER_ERROR_IF(not std::isfinite(width) or width < 0 or not std::isfinite(shift),
                         "Invalid ECS pressure width or shift for species {}",
                         spec)
      Wimag[k, k]               = width;
      real_val(Ws[page][k, k]) += weight * shift;
    }
    kernel(Wimag, bnd_qid, bnd, sort, spec, data->second, dipr, atm);
    sum_rule_residual[i] = closure_residual(Wimag, dipr);
    for (Size r = 0; r < n; ++r) {
      for (Size c = 0; c < n; ++c) { imag_val(Ws[page][r, c]) += weight * Wimag[r, c]; }
    }
    ++i;
  }
  for (Size i = 0; i < n; ++i) { Ws[joker, i, i] += bnd.lines[sort[i]].f0; }
}
ARTS_METHOD_ERROR_CATCH

void calculate(PropmatVectorView pm_,
               PropmatMatrixView,
               ComputeData&                    com_data,
               const ConstVectorView           f_grid_,
               const Range&                    f_range,
               const Jacobian::Targets&        jac_targets,
               const QuantumIdentifier&        bnd_qid,
               const band_data&                bnd,
               const LinemixingSpeciesEcsData& rovib_data,
               const AtmPoint&                 atm,
               const ZeemanPolarization        pol,
               const bool                      no_negative_absorption) try {
  if (pol != ZeemanPolarization::no) {
    ARTS_USER_ERROR_IF(stdr::any_of(
                           bnd, [](auto& zee) { return zee.on; }, &line::z),
                       "Zeeman effect and ECS in combination is not yet possible.")
    return;
  }

  PropmatVectorView     pm     = pm_[f_range];
  const ConstVectorView f_grid = f_grid_[f_range];

  ARTS_USER_ERROR_IF(jac_targets.target_count() > 0, "No Jacobian support.")

  if (bnd.size() == 0) return;

  com_data.adapt_single(bnd_qid, bnd, rovib_data, atm);

  com_data.core_calc(f_grid);

  for (Size i = 0; i < f_grid.size(); ++i) {
    const auto F =
        Constant::inv_sqrt_pi * atm[bnd_qid.isot.spec] * atm[bnd_qid.isot] * com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
  }
}
ARTS_METHOD_ERROR_CATCH

void equivalent_values(ComplexTensor3View              eqv_str,
                       ComplexTensor3View              eqv_val,
                       ComputeData&                    com_data,
                       const QuantumIdentifier&        bnd_qid,
                       const band_data&                bnd,
                       const LinemixingSpeciesEcsData& rovib_data,
                       const AtmPoint&                 atm,
                       const Vector&                   T) try {
  const auto k = eqv_str.npages();
  const auto m = eqv_str.ncols();

  ARTS_USER_ERROR_IF(eqv_str.shape() != eqv_val.shape(), "eqv_str and eqv_val must have the same shape.")
  ARTS_USER_ERROR_IF(T.size() != static_cast<Size>(k), "T must have the same size as eqv_str pages.")
  ARTS_USER_ERROR_IF(bnd.size() != static_cast<Size>(m), "bnd must have the same size as eqv_str cols.")

  if (bnd.size() == 0) return;
  ARTS_USER_ERROR_IF(eqv_str.nrows() != static_cast<Index>(bnd.front().ls.single_models.size()),
                     "eqv_str rows must match the number of ECS broadening species")

  com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm, false);

  std::string err{};
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(com_data)
  for (Index i = 0; i < k; ++i) {
    try {
      AtmPoint atm_copy    = atm;
      atm_copy.temperature = T[i];
      com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm_copy, true);
      com_data.core_calc_eqv();
      eqv_str[i] = com_data.eqv_strs;
      eqv_val[i] = com_data.eqv_vals;
    } catch (std::exception& e) {
#pragma omp critical
      err += std::format("{}\n", e.what());
    }
  }

  if (not err.empty()) throw std::runtime_error(err);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs
