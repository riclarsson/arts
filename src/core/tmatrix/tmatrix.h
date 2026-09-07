#pragma once

#include <rtepack.h>

namespace tmatrix {
//! Whether the optional, original Fortran backend is built.
bool available();
//! Whether the selected backend uses extended-precision internals.
bool extended_precision();

//! Owning results; no references to the Fortran COMMON blocks escape a call.
struct FixedResult {
  Index                           order{};
  Numeric                         scattering{}, extinction{};
  matpack::cdata_t<Complex, 2, 2> amplitude{};
  Muelmat                         phase{0.0};
};

//! Angles in degrees. Radius and wavelength must use the same length unit.
//! Amplitude has that length unit; phase and cross sections its square.
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
                  Numeric accuracy     = 0.001,
                  Numeric radius_ratio = 1,
                  int     shape        = -1);

struct RandomResult {
  Numeric effective_radius{}, effective_variance{}, extinction{}, scattering{}, albedo{}, asymmetry{};
  //! Rows at equally spaced scattering angles 0..180 degrees.
  //! Mueller matrices normalized as in Mishchenko (integral F11 = 4*pi).
  MuelmatVector phase;
};

//! One size distribution (NPNAX=1). Distribution parameters follow TMD.
//! Default is an effectively monodisperse, equal-volume-radius particle.
RandomResult random(Numeric radius,
                    Numeric wavelength,
                    Numeric aspect_ratio,
                    Numeric refractive_real,
                    Numeric refractive_imag,
                    int     angles             = 19,
                    Numeric accuracy           = 0.001,
                    Numeric radius_ratio       = 1,
                    int     shape              = -1,
                    int     distribution       = 4,
                    Numeric b                  = 0.1,
                    Numeric gamma              = 1,
                    int     size_quadrature    = -1,
                    int     surface_quadrature = 2,
                    Numeric lower_radius_ratio = 0.9999999,
                    Numeric upper_radius_ratio = 1.0000001);
}  // namespace tmatrix
