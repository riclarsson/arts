#include "tmatrix.h"

#include <cmath>
#include <mutex>
#include <stdexcept>
#include <string>

#ifdef ARTS_HAS_TMATRIX
#include <tmatrix_fortran.h>
#define tmatrix_ TMATRIX_FC_GLOBAL(tmatrix, TMATRIX)
#define ampl_    TMATRIX_FC_GLOBAL(ampl, AMPL)
#define tmd_     TMATRIX_FC_GLOBAL(tmd, TMD)
extern "C" {
void tmatrix_(double&,
              double&,
              int&,
              double&,
              double&,
              double&,
              double&,
              double&,
              int&,
              int&,
              double&,
              double&,
              char*,
              std::size_t);
void ampl_(int&, double&, double&, double&, double&, double&, double&, double&, Complex&, Complex&, Complex&, Complex&);
void tmd_(double&,
          int&,
          double&,
          int&,
          double&,
          double&,
          int&,
          double&,
          int&,
          double&,
          double&,
          double&,
          double&,
          int&,
          int&,
          double&,
          double&,
          int&,
          double&,
          double&,
          double&,
          double&,
          double&,
          double&,
          double*,
          double*,
          double*,
          double*,
          double*,
          double*,
          char*,
          std::size_t);
}
#endif

namespace tmatrix {
namespace {
#ifdef ARTS_HAS_TMATRIX
std::mutex solver_mutex;
#endif
void require(bool valid, const char* message) {
  if (!valid) throw std::invalid_argument(message);
}
void positive(double x, const char* name) { require(std::isfinite(x) && x > 0, name); }
void check(double radius,
           double wavelength,
           double aspect,
           double real,
           double imag,
           double accuracy,
           double ratio,
           int    shape) {
  if (!available()) throw std::runtime_error("T-matrix requires ENABLE_TMATRIX=ON");
  positive(radius, "radius must be finite and positive");
  positive(wavelength, "wavelength must be finite and positive");
  positive(aspect, "aspect_ratio must be finite and positive");
  positive(real, "refractive_real must be finite and positive");
  require(std::isfinite(imag) && imag >= 0, "refractive_imag must be finite and nonnegative");
  positive(accuracy, "accuracy must be finite and positive");
  positive(ratio, "radius_ratio must be finite and positive");
  require(shape == -1 || shape == -2, "shape must be -1 (spheroid) or -2 (cylinder)");
}
#ifdef ARTS_HAS_TMATRIX
void error(const char* msg) {
  std::string s(msg, 100);
  const auto  end = s.find_last_not_of(' ');
  if (end != std::string::npos) throw std::runtime_error("T-matrix: " + s.substr(0, end + 1));
}
void ampmat_to_phamat(Muelmat& z, const Complex& s11, const Complex& s12, const Complex& s21, const Complex& s22) {
  z[0, 0] = 0.5 * (s11 * conj(s11) + s12 * conj(s12) + s21 * conj(s21) + s22 * conj(s22)).real();
  z[0, 1] = 0.5 * (s11 * conj(s11) - s12 * conj(s12) + s21 * conj(s21) - s22 * conj(s22)).real();
  z[0, 2] = (-s11 * conj(s12) - s22 * conj(s21)).real();
  z[0, 3] = (Complex(0., 1.) * (s11 * conj(s12) - s22 * conj(s21))).real();

  z[1, 0] = 0.5 * (s11 * conj(s11) + s12 * conj(s12) - s21 * conj(s21) - s22 * conj(s22)).real();
  z[1, 1] = 0.5 * (s11 * conj(s11) - s12 * conj(s12) - s21 * conj(s21) + s22 * conj(s22)).real();
  z[1, 2] = (-s11 * conj(s12) + s22 * conj(s21)).real();
  z[1, 3] = (Complex(0., 1.) * (s11 * conj(s12) + s22 * conj(s21))).real();

  z[2, 0] = (-s11 * conj(s21) - s22 * conj(s12)).real();
  z[2, 1] = (-s11 * conj(s21) + s22 * conj(s12)).real();
  z[2, 2] = (s11 * conj(s22) + s12 * conj(s21)).real();
  z[2, 3] = (Complex(0., -1.) * (s11 * conj(s22) + s21 * conj(s12))).real();

  z[3, 0] = (Complex(0., 1.) * (s21 * conj(s11) + s22 * conj(s12))).real();
  z[3, 1] = (Complex(0., 1.) * (s21 * conj(s11) - s22 * conj(s12))).real();
  z[3, 2] = (Complex(0., -1.) * (s22 * conj(s11) - s12 * conj(s21))).real();
  z[3, 3] = (s22 * conj(s11) - s12 * conj(s21)).real();
}

#endif
}  // namespace

bool available() {
#ifdef ARTS_HAS_TMATRIX
  return true;
#else
  return false;
#endif
}

bool extended_precision() {
#ifdef ARTS_TMATRIX_QUAD
  return true;
#else
  return false;
#endif
}

FixedResult fixed(Numeric radius,
                  Numeric wavelength,
                  Numeric aspect_ratio,
                  Numeric refractive_real,
                  Numeric refractive_imag,
                  Numeric theta_incident,
                  Numeric theta_scattered,
                  Numeric phi_incident,
                  Numeric phi_scattered,
                  Numeric alpha,
                  Numeric beta,
                  Numeric accuracy,
                  Numeric radius_ratio,
                  int     shape) {
  matpack::cdata_t<Numeric, 1, 6> geometries;
  FixedResult                     result;
  geometries[0, 0] = theta_incident;
  geometries[0, 1] = theta_scattered;
  geometries[0, 2] = phi_incident;
  geometries[0, 3] = phi_scattered;
  geometries[0, 4] = alpha;
  geometries[0, 5] = beta;
  fixed_batch(std::span{&result, 1},
              radius,
              wavelength,
              aspect_ratio,
              refractive_real,
              refractive_imag,
              geometries,
              accuracy,
              radius_ratio,
              shape);
  return result;
}

void fixed_batch(std::span<FixedResult> results,
                 Numeric                radius,
                 Numeric                wavelength,
                 Numeric                aspect_ratio,
                 Numeric                refractive_real,
                 Numeric                refractive_imag,
                 ConstMatrixView        geometries,
                 Numeric                accuracy,
                 Numeric                radius_ratio,
                 int                    shape) {
  check(radius, wavelength, aspect_ratio, refractive_real, refractive_imag, accuracy, radius_ratio, shape);
  require(geometries.ncols() == 6, "geometries must have six columns");
  require(results.size() == static_cast<Size>(geometries.nrows()), "results must have one element per geometry row");
  for (Index i = 0; i < geometries.nrows(); ++i) {
    for (Index j : {0, 1, 5})
      require(std::isfinite(geometries[i, j]) && geometries[i, j] >= 0 && geometries[i, j] <= 180,
              "zenith angles and beta must be in [0, 180] degrees");
    for (Index j : {2, 3, 4})
      require(std::isfinite(geometries[i, j]) && geometries[i, j] >= 0 && geometries[i, j] <= 360,
              "azimuth angles and alpha must be in [0, 360] degrees");
  }
  if (geometries.nrows() == 0) return;
  // Equal-volume and equal-area radii coincide for a sphere. Avoid the
  // removable 0/0 in the legacy spheroid surface-area conversion.
  if (shape == -1 && aspect_ratio == 1) radius_ratio = 1;
#ifdef ARTS_HAS_TMATRIX
  FixedResult out;
  // AMPL reads the T-matrix from COMMON: keep the lock until all results
  // have been copied, including when alternating fixed/random calculations.
  std::lock_guard lock(solver_mutex);
  int             quiet = 1, order = 0;
  char            msg[100];
  tmatrix_(radius_ratio,
           radius,
           shape,
           wavelength,
           aspect_ratio,
           refractive_real,
           refractive_imag,
           accuracy,
           quiet,
           order,
           out.scattering,
           out.extinction,
           msg,
           sizeof(msg));
  error(msg);
  out.order = order;
  for (Index i = 0; i < geometries.nrows(); ++i) {
    auto theta_incident  = geometries[i, 0];
    auto theta_scattered = geometries[i, 1];
    auto phi_incident    = geometries[i, 2];
    auto phi_scattered   = geometries[i, 3];
    auto alpha           = geometries[i, 4];
    auto beta            = geometries[i, 5];
    ampl_(order,
          wavelength,
          theta_incident,
          theta_scattered,
          phi_incident,
          phi_scattered,
          alpha,
          beta,
          out.amplitude[0, 0],
          out.amplitude[0, 1],
          out.amplitude[1, 0],
          out.amplitude[1, 1]);
    ampmat_to_phamat(out.phase, out.amplitude[0, 0], out.amplitude[0, 1], out.amplitude[1, 0], out.amplitude[1, 1]);
    results[i] = out;
  }
#endif
}

std::vector<FixedResult> fixed_batch(Numeric         radius,
                                     Numeric         wavelength,
                                     Numeric         aspect_ratio,
                                     Numeric         refractive_real,
                                     Numeric         refractive_imag,
                                     ConstMatrixView geometries,
                                     Numeric         accuracy,
                                     Numeric         radius_ratio,
                                     int             shape) {
  std::vector<FixedResult> results(geometries.nrows());
  fixed_batch(results,
              radius,
              wavelength,
              aspect_ratio,
              refractive_real,
              refractive_imag,
              geometries,
              accuracy,
              radius_ratio,
              shape);
  return results;
}

RandomResult random(Numeric radius,
                    Numeric wavelength,
                    Numeric aspect_ratio,
                    Numeric refractive_real,
                    Numeric refractive_imag,
                    int     angles,
                    Numeric accuracy,
                    Numeric radius_ratio,
                    int     shape,
                    int     distribution,
                    Numeric b,
                    Numeric gamma,
                    int     size_quadrature,
                    int     surface_quadrature,
                    Numeric lower_radius_ratio,
                    Numeric upper_radius_ratio) {
  check(radius, wavelength, aspect_ratio, refractive_real, refractive_imag, accuracy, radius_ratio, shape);
  require(angles >= 2, "angles must be at least 2");
  require(distribution >= 1 && distribution <= 5, "distribution must be in [1, 5]");
  if (distribution == 5)
    require(std::isfinite(b), "b must be finite");
  else
    positive(b, "b must be finite and positive");
  positive(gamma, "gamma must be finite and positive");
  require(size_quadrature >= -1 && size_quadrature <= 998, "size_quadrature must be in [-1, 998]");
  if (distribution == 5) require(size_quadrature >= 1, "modified power law requires size_quadrature >= 1");
  require(surface_quadrature >= 1 && surface_quadrature <= 900, "surface_quadrature must be in [1, 900]");
  positive(lower_radius_ratio, "lower_radius_ratio must be finite and positive");
  require(std::isfinite(upper_radius_ratio) && upper_radius_ratio > lower_radius_ratio,
          "upper_radius_ratio must exceed lower_radius_ratio");
  if (shape == -1 && aspect_ratio == 1) radius_ratio = 1;
  // These distributions parameterize radius by a length. Normalizing it
  // avoids the absolute 1e-5 root bracket in the legacy POWER routine.
  const Numeric length_scale  = (distribution == 3 || distribution == 4) ? radius : 1;
  radius                     /= length_scale;
  wavelength                 /= length_scale;
  RandomResult out;
  out.phase.resize(angles);
#ifdef ARTS_HAS_TMATRIX
  std::lock_guard lock(solver_mutex);
  int             count = 1, quiet = 1;
  char            msg[100];
  // Fortran outputs six contiguous vectors; transpose to ARTS row-major rows.
  Matrix phase(6, angles);
  tmd_(radius_ratio,
       distribution,
       radius,
       count,
       b,
       gamma,
       size_quadrature,
       aspect_ratio,
       shape,
       wavelength,
       refractive_real,
       refractive_imag,
       accuracy,
       angles,
       surface_quadrature,
       lower_radius_ratio,
       upper_radius_ratio,
       quiet,
       out.effective_radius,
       out.effective_variance,
       out.extinction,
       out.scattering,
       out.albedo,
       out.asymmetry,
       &phase[0, 0],
       &phase[1, 0],
       &phase[2, 0],
       &phase[3, 0],
       &phase[4, 0],
       &phase[5, 0],
       msg,
       sizeof(msg));
  error(msg);
  out.effective_radius *= length_scale;
  out.extinction       *= length_scale * length_scale;
  out.scattering       *= length_scale * length_scale;
  for (int i = 0; i < angles; ++i)
    out.phase[i] = Muelmat{phase[0, i],
                           phase[4, i],
                           0,
                           0,
                           phase[4, i],
                           phase[1, i],
                           0,
                           0,
                           0,
                           0,
                           phase[2, i],
                           phase[5, i],
                           0,
                           0,
                           -phase[5, i],
                           phase[3, i]};
#endif
  return out;
}
}  // namespace tmatrix
