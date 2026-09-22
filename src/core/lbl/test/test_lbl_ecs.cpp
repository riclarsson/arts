#include <arts_constants.h>
#include <atm.h>
#include <jacobian.h>
#include <lbl_lineshape_voigt_ecs.h>
#include <lbl_lineshape_voigt_ecs_makarov.h>
#include <lbl_lineshape_voigt_lte.h>
#include <physics_funcs.h>
#include <wigner_functions.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
void require(bool ok, std::string_view message) {
  if (not ok) throw std::runtime_error(std::string(message));
}
void near(Complex actual, Complex expected, Numeric tolerance, std::string_view message) {
  require(std::abs(actual - expected) <= tolerance * std::max(std::abs(expected), Numeric{1e-100}), message);
}
template <class F> void throws(F&& f, std::string_view text) {
  try {
    f();
  } catch (const std::exception& e) {
    require(std::string_view(e.what()).contains(text), e.what());
    return;
  }
  throw std::runtime_error(std::format("Expected an error containing '{}'", text));
}

QuantumIdentifier co2_id() {
  QuantumIdentifier id{"CO2-626"_isot};
  id.state[QuantumNumberType::l2] = {.upper = Rational{0}, .lower = Rational{0}};
  return id;
}

lbl::band_data co2_band() {
  lbl::band_data band;
  band.lineshape = LineByLineLineshape::VP_ECS_HARTMANN;
  lbl::line ln;
  ln.f0                       = 1e11;
  ln.a                        = 1e-5;
  ln.gu                       = 3;
  ln.gl                       = 1;
  ln.ls.T0                    = 296;
  ln.qn[QuantumNumberType::J] = {.upper = Rational{1}, .lower = Rational{0}};
  ln.ls.single_models[SpeciesEnum::Nitrogen].data[LineShapeModelVariable::G0] =
      lbl::temperature::data{LineShapeModelType::T0, {1.0}};
  band.lines.push_back(ln);
  return band;
}

LinemixingSpeciesEcsData collision_data() {
  LinemixingSingleEcsData rate;
  rate.scaling              = {LineShapeModelType::T0, {1.0}};
  rate.beta                 = {LineShapeModelType::T0, {0.01}};
  rate.lambda               = {LineShapeModelType::T0, {0.6}};
  rate.collisional_distance = {LineShapeModelType::T0, {2e-10}};
  return {{SpeciesEnum::Nitrogen, rate}, {SpeciesEnum::Oxygen, rate}};
}

void isolated_line() {
  const auto id    = co2_id();
  const auto band  = co2_band();
  const auto rates = collision_data();
  for (Numeric pressure : {1.0, 1e5, 1e7}) {
    AtmPoint atm;
    atm.pressure                    = pressure;
    atm.temperature                 = 296;
    atm[SpeciesEnum::Nitrogen]      = 1;
    atm[SpeciesEnum::CarbonDioxide] = 4e-4;
    atm[id.isot]                    = 1;
    const Numeric D =
        band.front().f0 * std::sqrt(Constant::doppler_broadening_const_squared * atm.temperature / id.isot.mass);
    Vector                       grid{band.front().f0 - 2 * D, band.front().f0, band.front().f0 + D};
    lbl::voigt::ecs::ComputeData data(grid, atm);
    PropmatVector                pm(grid.size());
    pm = 0;
    PropmatMatrix     dpm(0, grid.size());
    Jacobian::Targets targets;
    lbl::voigt::ecs::calculate(
        pm, dpm, data, grid, Range(0, grid.size()), targets, id, band, rates, atm, ZeemanPolarization::no, false);
    lbl::voigt::lte::single_shape isolated(id.isot, band.front(), atm, ZeemanPolarization::no, 0);
    for (Size i = 0; i < grid.size(); ++i) {
      const Numeric factor = data.scl[i] * Constant::c * Constant::c / (8 * Constant::pi);
      near(pm[i].A(), factor * isolated(grid[i]).real(), 2e-12, "Isolated ECS differs from ordinary Voigt");
    }
    require(data.sum_rule_residual[0] == 1, "Single-line closure residual should expose the missing couplings");
  }
}

void eigen_resolvent() {
  AtmPoint atm;
  atm.pressure    = 1e5;
  atm.temperature = 296;
  Vector                       grid{1e9};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  data.pop  = Vector{2.0, 0.5};
  data.dip  = Vector{1.0, -0.8};
  data.vmrs = Vector{1.0};
  data.Ws.resize(1, 2, 2);
  data.Ws[0][0, 0] = Complex(1e9, 2);
  data.Ws[0][1, 1] = Complex(1e9 + 5, 3);
  data.Ws[0][0, 1] = Complex(0, -0.2);
  data.Ws[0][1, 0] = Complex(0, -0.1);
  const ComplexMatrix original{data.Ws[0]};
  data.core_calc_eqv();
  const ComplexMatrix strengths{data.eqv_strs};
  const ComplexMatrix positions{data.eqv_vals};
  data.core_calc_eqv();
  Complex total = 0;
  for (Index i = 0; i < 2; ++i) {
    near(data.eqv_strs[0, i], strengths[0, i], 1e-14, "Equivalent strengths changed on repeated calculation");
    near(data.eqv_vals[0, i], positions[0, i], 1e-14, "Equivalent positions changed on repeated calculation");
    total += data.eqv_strs[0, i];
    for (Index j = 0; j < 2; ++j) require(data.Ws[0][i, j] == original[i, j], "ECS input matrix mutated");
  }
  near(total, 2.0 + 0.5 * 0.8 * 0.8, 2e-14, "Equivalent strengths violate residue-sum invariant");
  require(data.eigenvector_rcond[0] > 0.1, "Unexpectedly ill-conditioned test eigenvectors");
  for (Numeric frequency : {1e9 - 4, 1e9 + 2, 1e9 + 9}) {
    ComplexMatrix A(2, 2);
    ComplexVector rhs(2), solution(2);
    for (Index i = 0; i < 2; ++i) {
      rhs[i] = data.pop[i] * data.dip[i];
      for (Index j = 0; j < 2; ++j) A[i, j] = (i == j ? Complex(frequency) : Complex(0)) - original[j, i];
    }
    solve(solution, A, rhs);
    Complex direct = 0, equivalent = 0;
    for (Index i = 0; i < 2; ++i) {
      direct     += data.dip[i] * solution[i];
      equivalent += data.eqv_strs[0, i] / (frequency - data.eqv_vals[0, i]);
    }
    near(equivalent, direct, 5e-8, "Equivalent-line resolvent differs from direct solve");
  }
  data.Ws[0][0, 0] = Complex(1e9, 2);
  data.Ws[0][1, 1] = Complex(1e9, 2);
  data.Ws[0][0, 1] = Complex(0, 1);
  data.Ws[0][1, 0] = 0;
  throws([&] { data.core_calc_eqv(); }, "condition");
  data.Ws[0][0, 1] = 0;
  data.Ws[0][0, 0] = Complex(1e9, -1);
  throws([&] { data.core_calc_eqv(); }, "negative damping");
}

void makarov_and_presorting() {
  QuantumIdentifier id{"O2-66"_isot};
  id.state[QuantumNumberType::S]        = {.upper = Rational{1}, .lower = Rational{1}};
  auto band                             = co2_band();
  band.lineshape                        = LineByLineLineshape::VP_ECS_MAKAROV;
  band.front().qn[QuantumNumberType::N] = {.upper = Rational{1}, .lower = Rational{1}};
  AtmPoint atm;
  atm.pressure               = 1e5;
  atm.temperature            = 296;
  atm[SpeciesEnum::Nitrogen] = 1;
  Vector                       grid{band.front().f0};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  data.adapt_single(id, band, collision_data(), atm);
  data.core_calc(grid);
  require(data.shape[0].real() > 0, "Supported Makarov single line has no absorption");
  atm.specs.clear();
  id.state[QuantumNumberType::S].upper = Rational{0};
  throws([&] { data.adapt_single(id, band, {}, atm); }, "S=1");
  id.state[QuantumNumberType::S].upper        = Rational{1};
  band.front().qn[QuantumNumberType::N].lower = Rational{3};
  throws([&] { data.adapt_single(id, band, {}, atm); }, "N");

  band = co2_band();
  for (Index i = 1; i < 3; ++i) {
    auto ln                      = band.front();
    ln.f0                       += Numeric(i) * 1e7;
    ln.a                        *= Numeric(i + 1);
    ln.qn[QuantumNumberType::J]  = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i}};
    band.lines.push_back(ln);
  }
  const auto co2 = co2_id();
  Vector     expected_pop(3), expected_dip(3), expected_dipr(3);
  for (Index i = 0; i < 3; ++i) {
    lbl::band_data single;
    single.lineshape = band.lineshape;
    single.lines.push_back(band.lines[i]);
    data.adapt_single(co2, single, {}, atm);
    expected_pop[i]  = data.pop[0];
    expected_dip[i]  = data.dip[0];
    expected_dipr[i] = data.dipr[0];
  }
  data.sort = ArrayOfIndex{1, 2, 0};
  data.adapt_single(co2, band, {}, atm, true);
  for (Index i = 0; i < 3; ++i) {
    const Index j = data.sort[i];
    near(data.pop[i], expected_pop[j], 1e-14, "Presorted population mismatch");
    near(data.dip[i], expected_dip[j], 1e-14, "Presorted dipole mismatch");
    near(data.dipr[i], expected_dipr[j], 1e-14, "Presorted reduced dipole mismatch");
    near(data.Ws[0][i, i].real(), band.lines[j].f0, 1e-14, "Presorted line position mismatch");
  }
}

void validation_and_mixtures() {
  auto       band  = co2_band();
  const auto id    = co2_id();
  const auto rates = collision_data();
  AtmPoint   atm;
  atm.pressure               = 1e5;
  atm.temperature            = 296;
  atm[SpeciesEnum::Nitrogen] = 0.8;
  atm[SpeciesEnum::Oxygen]   = 0.2;
  Vector                       grid{1e11};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  auto                         oxygen                = band.front().ls.single_models.at(SpeciesEnum::Nitrogen);
  oxygen.data[LineShapeModelVariable::G0]            = {LineShapeModelType::T0, {3.0}};
  oxygen.data[LineShapeModelVariable::D0]            = {LineShapeModelType::T0, {0.1}};
  band.front().ls.single_models[SpeciesEnum::Oxygen] = oxygen;
  data.adapt_single(id, band, rates, atm);
  const Complex mixed = data.Ws[0][0, 0];
  data.adapt_multi(id, band, rates, atm);
  Complex weighted = 0;
  for (Size i = 0; i < data.vmrs.size(); ++i) weighted += data.vmrs[i] * data.Ws[i][0, 0];
  near(mixed, weighted, 1e-14, "Single ECS mixture differs from sum of partner matrices");
  near(mixed.imag(), 1.4e5, 1e-14, "Wrong mixture broadening");

  // Unordered maps with the same keys need not have the same iteration order.
  auto next                      = band.front();
  next.f0                       += 1e7;
  next.qn[QuantumNumberType::J]  = {.upper = Rational{3}, .lower = Rational{2}};
  next.ls.single_models.clear();
  next.ls.single_models[SpeciesEnum::Oxygen]   = oxygen;
  next.ls.single_models[SpeciesEnum::Nitrogen] = band.front().ls.single_models.at(SpeciesEnum::Nitrogen);
  band.lines.push_back(next);
  data.adapt_single(id, band, rates, atm);
  band.lines.back().ls.single_models.erase(SpeciesEnum::Oxygen);
  throws([&] { data.adapt_single(id, band, rates, atm); }, "same broadening species");
  throws([&] { data.adapt_multi(id, band, rates, atm); }, "same broadening species");

  band = co2_band();
  band.lines.push_back(band.front());
  throws([&] { data.adapt_single(id, band, rates, atm); }, "same rotational pair");
  band = co2_band();
  atm.specs.clear();
  auto unsupported_id = id;
  unsupported_id.isot = "CO2-628"_isot;
  throws([&] { data.adapt_single(unsupported_id, band, {}, atm); }, "CO2-626");
  data.adapt_single(id, band, {}, atm);
  require(data.Ws[0][0, 0].imag() == 0, "Absent perturbers should produce zero width");
  data.core_calc(grid);
  require(std::isfinite(data.shape[0].real()), "Absent perturbers produced non-finite spectrum");
  atm[SpeciesEnum::Nitrogen] = -0.1;
  throws([&] { data.adapt_single(id, band, rates, atm); }, "VMR");

  ComplexTensor3 strength(1, 1, 1), position(2, 1, 1);
  Vector         temperatures{296};
  throws([&] { lbl::voigt::ecs::equivalent_values(strength, position, data, id, band, rates, atm, temperatures); },
         "same shape");
}

void explicit_partner_mixture() {
  const auto id   = co2_id();
  auto       band = co2_band();
  for (Index i = 1; i < 3; ++i) {
    auto ln                      = band.front();
    ln.f0                       += Numeric(i) * 1e7;
    ln.a                        *= Numeric(i + 1);
    ln.qn[QuantumNumberType::J]  = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i}};
    band.lines.push_back(ln);
  }
  for (Size i = 0; i < band.size(); ++i) {
    auto& models                            = band.lines[i].ls.single_models;
    auto  oxygen                            = models.at(SpeciesEnum::Nitrogen);
    oxygen.data[LineShapeModelVariable::G0] = {LineShapeModelType::T0, {2.0 + Numeric(i)}};
    oxygen.data[LineShapeModelVariable::D0] = {LineShapeModelType::T0, {0.1 * Numeric(i + 1)}};
    models[SpeciesEnum::Oxygen]             = oxygen;
  }
  auto rates                                         = collision_data();
  rates.at(SpeciesEnum::Oxygen).beta                 = {LineShapeModelType::T0, {0.05}};
  rates.at(SpeciesEnum::Oxygen).lambda               = {LineShapeModelType::T0, {0.8}};
  rates.at(SpeciesEnum::Oxygen).collisional_distance = {LineShapeModelType::T0, {3e-10}};
  AtmPoint atm;
  atm.pressure               = 1e5;
  atm[SpeciesEnum::Nitrogen] = 0.79;
  atm[SpeciesEnum::Oxygen]   = 0.21;
  Vector grid{band.front().f0};
  for (Numeric temperature : {220.0, 296.0, 330.0}) {
    atm.temperature = temperature;
    lbl::voigt::ecs::ComputeData data(grid, atm);
    data.adapt_single(id, band, rates, atm);
    const ComplexMatrix mixed{data.Ws[0]};
    data.adapt_multi(id, band, rates, atm);
    bool has_coupling = false;
    for (Size r = 0; r < band.size(); ++r) {
      for (Size c = 0; c < band.size(); ++c) {
        Complex expected = 0;
        for (Size p = 0; p < data.vmrs.size(); ++p) expected += data.vmrs[p] * data.Ws[p][r, c];
        near(mixed[r, c], expected, 2e-14, "ECS mixture must sum the full collision-partner matrices");
        has_coupling |= r != c and std::abs(mixed[r, c]) > 0;
      }
    }
    require(has_coupling, "Mixture fixture must exercise off-diagonal ECS couplings");
  }
}
}  // namespace

int main() try {
  const auto rates = collision_data();
  throws([&] { (void)rates.at(SpeciesEnum::Nitrogen).Q(Rational{0}, 296, 296, 0); }, "Q(0)");
  WignerInformation wigner(100, 0, true, true);
  isolated_line();
  eigen_resolvent();
  makarov_and_presorting();
  validation_and_mixtures();
  explicit_partner_mixture();
  wigner.unload();
  std::cout << "ECS numerical regression tests passed\n";
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
