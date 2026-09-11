#include <disort.h>

#include <cmath>
#include <stdexcept>

int main() {
  constexpr Index retained_moments = 2;
  constexpr Index full_moments     = 4;

  DisortSettings settings;
  settings.resize(2, retained_moments, 1, AscendingGrid{1.0}, DescendingGrid{1000.0, 0.0});
  settings.legendre_coefficients.resize(1, 1, full_moments);
  settings.legendre_coefficients[0, 0, joker] = Vector{1.0, 0.9, 0.8, 0.72};

  const auto scaling                = disort::delta_m_plus(settings.legendre_coefficients[0], retained_moments);
  settings.fractional_scattering[0] = scaling.fraction;
  settings.delta_m_peak_moments[0]  = scaling.moments;

  const auto solver = settings.init();
  if (solver.delta_m_peak_moments().ncols() != retained_moments)
    throw std::runtime_error("Delta-M-plus retained moments were not assigned to NLeg");
  if (solver.all_legendre_coeffs().ncols() != full_moments)
    throw std::runtime_error("Delta-M-plus full moments were not assigned to NLeg_all");

  // DFDT must differentiate total net flux in physical optical depth even
  // when the beam is delta-M scaled. Check all three scalar flux interfaces.
  settings.resize(4, 4, 1, AscendingGrid{1.0}, DescendingGrid{1000.0, 0.0});
  settings.optical_thicknesses                = 1.0;
  settings.single_scattering_albedo           = 0.6;
  settings.fractional_scattering              = 0.2;
  settings.legendre_coefficients[0, 0, joker] = Vector{1.0, 0.7, 0.5, 0.3};
  settings.solar_source                       = 2.0;
  settings.solar_zenith_angle                 = 60.0;
  settings.source_polynomial.resize(1, 1, 1);
  settings.source_polynomial = 1.0;
  auto dis                   = settings.init();
  settings.set(dis, 0);
  disort::flux_data scratch;
  const auto        net = [&](Numeric tau) {
    const auto f = dis.flux(scratch, tau);
    return f.up - f.down_diffuse - f.down_direct;
  };
  constexpr Numeric h                   = 1e-5;
  const Numeric     interior_derivative = (net(0.4 + h) - net(0.4 - h)) / (2 * h);
  const Numeric     boundary_derivative = (3 * net(1.0) - 4 * net(1.0 - h) + net(1.0 - 2 * h)) / (2 * h);
  const auto        check               = [](Numeric actual, Numeric expected) {
    if (not std::isfinite(actual) or std::abs(actual - expected) > 1e-8)
      throw std::runtime_error("Delta-M DFDT does not differentiate net flux in physical optical depth");
  };
  check(dis.flux(scratch, 0.4).dfdt, interior_derivative);
  check(dis.flux(scratch, 1.0).dfdt, boundary_derivative);
  Vector up(2), down(2), direct(2), dfdt(2);
  dis.ungridded_flux(up, down, direct, dfdt, AscendingGrid{0.4, 1.0});
  check(dfdt[0], interior_derivative);
  check(dfdt[1], boundary_derivative);
  up.resize(1);
  down.resize(1);
  direct.resize(1);
  dfdt.resize(1);
  dis.gridded_flux(up, down, direct, dfdt);
  check(dfdt[0], boundary_derivative);
}
