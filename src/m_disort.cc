#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <format>
#include <numeric>
#include <vector>

#include "arts_constants.h"
#include "configtypes.h"
#include "debug.h"
#include "path_point.h"

////////////////////////////////////////////////////////////////////////
// Core Disort
////////////////////////////////////////////////////////////////////////

void disort_spectral_radiance_fieldCalc(Tensor4& disort_spectral_radiance_field,
                                        Vector& disort_quadrature_angles,
                                        Vector& disort_quadrature_weights,
                                        const DisortSettings& disort_settings,
                                        const Vector& phis) {
  ARTS_TIME_REPORT

  using Conversion::acosd;

  const Index nv    = disort_settings.nfreq;
  const Index np    = disort_settings.nlay;
  const Index nquad = disort_settings.quadrature_dimension;

  disort_spectral_radiance_field.resize(nv, np, phis.size(), nquad);

  disort::main_data dis = disort_settings.init();
  Tensor3 ims(np, phis.size(), nquad / 2);
  Tensor3 tms(np, phis.size(), nquad);

  //! Supplementary outputs
  disort_quadrature_weights = dis.weights();
  disort_quadrature_angles.resize(nquad);
  std::transform(dis.mu().begin(),
                 dis.mu().end(),
                 disort_quadrature_angles.begin(),
                 [](const Numeric& mu) { return acosd(mu); });

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(dis, tms, ims)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_settings.set(dis, iv);
      dis.gridded_u_corr(disort_spectral_radiance_field[iv], tms, ims, phis);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(
      error.size(), "Error occurred in disort-spectral:\n{}", error);
}

void disort_spectral_flux_fieldCalc(Tensor3& disort_spectral_flux_field,
                                    const DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  const Index nv = disort_settings.nfreq;
  const Index np = disort_settings.nlay;

  disort_spectral_flux_field.resize(nv, 3, np);

  disort::main_data dis = disort_settings.init();

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_settings.set(dis, iv);

      dis.gridded_flux(disort_spectral_flux_field[iv, 0, joker],
                       disort_spectral_flux_field[iv, 1, joker],
                       disort_spectral_flux_field[iv, 2, joker]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n{}", error);
}

////////////////////////////////////////////////////////////////////////
// Integrate functions
////////////////////////////////////////////////////////////////////////

void spectral_radianceIntegrateDisort(
    StokvecVector& /*spectral_radiance*/,
    const Tensor4& /*disort_spectral_radiance_field*/,
    const Vector& /*disort_quadrature_angles*/,
    const Vector& /*disort_quadrature_weights*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}

void SpectralFluxDisort(Matrix& /*spectral_flux_field_up*/,
                        Matrix& /*spectral_flux_field_down*/,
                        const Tensor3& /*disort_spectral_flux_field*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}
