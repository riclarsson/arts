#include "disort.h"

#include <matpack.h>

#include <algorithm>
#include <cstring>
#include <ranges>
#include <vector>

#include "arts_constants.h"
#include "compare.h"
#include "legendre.h"
#include "lin_alg.h"

// #define TIMEIT

#ifdef TIMEIT
#include <artstime.h>

#include <iostream>
#define TIMEMACRO(x) x
#else
#define TIMEMACRO(x)
#endif

namespace disort {
void mathscr_v(auto&& um,
               mathscr_v_data& data,
               const Numeric tau,
               const ExhaustiveConstVectorView& source_poly_coeffs,
               const ExhaustiveConstMatrixView& G,
               const ExhaustiveConstVectorView& K,
               const ExhaustiveConstVectorView& inv_mu,
               const Index Ni0 = 0,
               const Numeric scl = 1.0,
               const Numeric add = 0.0) {
  const Index Ni = um.size();
  const Index Nk = G.ncols();
  const Index Nc = source_poly_coeffs.size();

  for (Index c = 0; c < Nc; c++) {
    data.cvec[c] = std::pow(tau, Nc - 1 - c);
  }

  data.src = 0.0;
  for (Index k = 0; k < Nk; k++) {
    for (Index i : std::ranges::iota_view(0, Nc)) {
      for (Index j = 0; j <= i; j++) {
        const Numeric fs =
            Legendre::tgamma_ratio(static_cast<Numeric>(Nc - j),
                                   static_cast<Numeric>(Nc - i)) *
            source_poly_coeffs[1 - j];
        data.src(k, i) += fs * std::pow(K[k], -(i - j + 1));
      }
    }
  }

  std::ranges::copy(inv_mu, data.k1.begin());
  std::copy(G.elem_begin(), G.elem_end(), data.G.elem_begin());
  solve_inplace(data.k1, data.G, data.solve_work);
  mult(data.k2, data.src, data.cvec);
  data.k1 *= data.k2;
  mult(um, G.slice(Ni0, Ni), data.k1, scl, add);
}

void main_data::solve_for_coefs() {
  const Index ln = NLayers - 1;
  auto RHS_middle = RHS.slice(N, n - NQuad).reshape_as(NQuad, NLayers - 1);

  TIMEMACRO(Numeric dtm{});
  TIMEMACRO(Numeric dtcollec{});
  TIMEMACRO(Numeric dtsolv{});
  TIMEMACRO(Numeric dtlhs{});
  TIMEMACRO(Numeric dtrhs{});
  TIMEMACRO(Numeric dtstart{});

  TIMEMACRO(Time total{});
  for (Index m = 0; m < NFourier; m++) {
    TIMEMACRO(Time tm{});
    const bool m_equals_0_bool = m == 0;
    const bool BDRF_bool = m < NBDRF;
    const auto G_collect_m = G_collect[m];
    const auto K_collect_m = K_collect[m];
    const auto B_collect_m = B_collect[m];

    if (BDRF_bool) {
      brdf_fourier_modes[m](
          mathscr_D_neg, mu_arr.slice(0, N), mu_arr.slice(N, N)),
          einsum<"ij", "", "ij", "j", "j">(
              R, 1 + m_equals_0_bool, mathscr_D_neg, mu_arr.slice(0, N), W);
      if (has_beam_source) {
        brdf_fourier_modes[m](mathscr_X_pos.reshape_as(N, 1),
                              mu_arr.slice(0, N),
                              ExhaustiveConstVectorView{-mu0});
        mathscr_X_pos *= mu0 * I0 / Constant::pi;
      }
    }

    const auto b_pos_m = b_pos[m];
    const auto b_neg_m = b_neg[m];

    TIMEMACRO(dtstart += TimeStep(Time{} - tm).count());

    // Fill RHS
    {
      TIMEMACRO(Time trhs{});
      if (has_source_poly and m_equals_0_bool) {
        mathscr_v(RHS.slice(0, N),
                  comp_data,
                  0.0,
                  source_poly_coeffs[0],
                  G_collect_m[0],
                  K_collect_m[0],
                  inv_mu_arr,
                  N,
                  -1.0);

        if (is_multilayer) {
          for (Index l = 0; l < ln; l++) {
            mathscr_v(RHS(Range(l + N, NQuad, ln)),
                      comp_data,
                      tau_arr[l],
                      source_poly_coeffs[l],
                      G_collect_m[l],
                      K_collect_m[l],
                      inv_mu_arr);

            mathscr_v(RHS(Range(l + N, NQuad, ln)),
                      comp_data,
                      tau_arr[l],
                      source_poly_coeffs[l + 1],
                      G_collect_m[l + 1],
                      K_collect_m[l + 1],
                      inv_mu_arr,
                      0,
                      -1.0,
                      1.0);
          }
        }

        mathscr_v(RHS.slice(n - N, N),
                  comp_data,
                  tau_arr.back(),
                  source_poly_coeffs[ln],
                  G_collect_m[ln],
                  K_collect_m[ln],
                  inv_mu_arr,
                  N);

        if (NBDRF > 0) {
          jvec = RHS.slice(n - N, N);
          mult(RHS.slice(n - N, N), R, jvec, 1.0, 1.0);
        }
      } else {
        RHS = 0.0;
      }

      if (has_beam_source) {
        if (BDRF_bool) {
          std::ranges::copy(mathscr_X_pos, BDRF_RHS_contribution.begin());
          mult(BDRF_RHS_contribution, R, B_collect_m[ln].slice(0, N), 1.0, 1.0);
        } else {
          BDRF_RHS_contribution = 0.0;
        }

        if (is_multilayer) {
          for (Index l = 0; l < ln; l++) {
            const Numeric scl = std::exp(-mu0 * scaled_tau_arr_with_0[l + 1]);
            for (Index j = 0; j < NQuad; j++) {
              RHS_middle(j, l) +=
                  (B_collect_m(l + 1, j) - B_collect_m(l, j)) * scl;
            }
          }
        }

        for (Index i = 0; i < N; i++) {
          RHS[i] += b_neg_m[i] - B_collect_m(0, N + i);
          RHS[n - N + i] +=
              b_pos_m[i] + (BDRF_RHS_contribution[i] - B_collect_m(ln, i)) *
                               std::exp(-scaled_tau_arr_with_0.back() / mu0);
        }
      } else {
        RHS.slice(0, N) += b_neg_m;
        RHS.slice(n - N, N) += b_pos_m;
      }

      TIMEMACRO(dtrhs += TimeStep(Time{} - trhs).count());
    }

    // Fill LHS
    {
      TIMEMACRO(Time dlhs{});
      for (Index i = 0; i < N; i++) {
        E_Lm1L[i] =
            std::exp(K_collect_m(ln, i) * (scaled_tau_arr_with_0[NLayers] -
                                           scaled_tau_arr_with_0[NLayers - 1]));
      }

      if (BDRF_bool) {
        mult(BDRF_LHS, R, G_collect_m[ln].slice(N, N));
      } else {
        BDRF_LHS = 0;
      }

      for (Index i = 0; i < N; i++) {
        for (Index j = 0; j < N; j++) {
          LHSB(i, j) = G_collect_m(0, i + N, j);
          LHSB(i, N + j) =
              G_collect_m(0, i + N, j + N) *
              std::exp(K_collect_m(0, j) * scaled_tau_arr_with_0[1]);
          LHSB(n - N + i, n - 2 * N + j) =
              (G_collect_m(ln, i, j) - BDRF_LHS(i, j)) * E_Lm1L[j];
          LHSB(n - N + i, n - N + j) =
              G_collect_m(ln, i, j + N) - BDRF_LHS(i, j + N);
        }
      }

      for (Index l = 0; l < ln; l++) {
        for (Index i = 0; i < N; i++) {
          E_lm1l[i] =
              std::exp(K_collect_m(l, i + N) * (scaled_tau_arr_with_0[l] -
                                                scaled_tau_arr_with_0[l + 1]));
          E_llp1[i] = std::exp(
              K_collect_m(l + 1, i + N) *
              (scaled_tau_arr_with_0[l + 1] - scaled_tau_arr_with_0[l + 2]));
        }

        for (Index i = 0; i < N; i++) {
          for (Index j = 0; j < N; j++) {
            LHSB(N + l * NQuad + i, l * NQuad + j) =
                G_collect_m(l, i, j) * E_lm1l[j];
            LHSB(2 * N + l * NQuad + i, l * NQuad + j) =
                G_collect_m(l, N + i, j) * E_lm1l[j];
            LHSB(N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j) =
                -G_collect_m(l + 1, i, N + j) * E_llp1[j];
            LHSB(2 * N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j) =
                -1 * G_collect_m(l + 1, N + i, N + j) * E_llp1[j];
          }
        }

        for (Index i = 0; i < NQuad; i++) {
          for (Index j = 0; j < N; j++) {
            LHSB(N + l * NQuad + i, l * NQuad + N + j) =
                G_collect_m(l, i, N + j);
            LHSB(N + l * NQuad + i, l * NQuad + 2 * N + j) =
                -G_collect_m(l + 1, i, j);
          }
        }
      }
      TIMEMACRO(dtlhs += TimeStep(Time{} - dlhs).count());
    }

    TIMEMACRO(Time solv{});
    LHSB.solve(RHS);
    TIMEMACRO(dtsolv += TimeStep(Time{} - solv).count());

    TIMEMACRO(Time collec{});
    einsum<"ijm", "ijm", "im">(
        GC_collect[m], G_collect_m, RHS.reshape_as(NLayers, NQuad));
    TIMEMACRO(dtcollec += TimeStep(Time{} - collec).count());

    TIMEMACRO(dtm += TimeStep(Time{} - tm).count());
  }

  TIMEMACRO(const auto ratio = [](auto a, auto b) {
    return std::round(1000 * a / b) / 10;
  });
  TIMEMACRO(std::cout << "Total solve-for-coeffs-time: " << dtm << '\n');
  TIMEMACRO(std::cout << "Ratio in START: " << ratio(dtstart, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in RHS:   " << ratio(dtrhs, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in LHS:   " << ratio(dtlhs, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in SOLVE: " << ratio(dtsolv, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in GCCOL: " << ratio(dtcollec, dtm) << '%'
                      << '\n');
}

Numeric poch(Numeric x, Numeric n) { return Legendre::tgamma_ratio(x + n, x); }

void main_data::diagonalize() {
  auto GmG = Gml.slice(0, N);

  TIMEMACRO(Numeric dtm{});
  TIMEMACRO(Numeric dtstart{});
  TIMEMACRO(Numeric dtspecial{});
  TIMEMACRO(Numeric dtlstart{});
  TIMEMACRO(Numeric dtbeam{});
  TIMEMACRO(Numeric dtD{});
  TIMEMACRO(Numeric dtsqr{});
  TIMEMACRO(Numeric dtdiag{});
  TIMEMACRO(Numeric dtcollectG{});

  for (Index m = 0; m < NFourier; m++) {
    TIMEMACRO(Time tm{});
    auto Km = K_collect[m];
    auto Gm = G_collect[m];
    auto Bm = B_collect[m];

    D_temp.resize(N, NLeg - m);
    X_temp.resize(NLeg - m);
    asso_leg_term_pos.resize(NLeg - m, N);
    asso_leg_term_neg.resize(NLeg - m, N);
    asso_leg_term_mu0.resize(NLeg - m);
    weighted_asso_Leg_coeffs_l.resize(NLeg - m);

    auto xtemp = X_temp.reshape_as(1, NLeg - m);

    const bool m_equals_0_bool = (m == 0);

    fac.resize(NLeg - m);
    for (Index i = m; i < NLeg; i++) {
      fac[i - m] =
          poch(static_cast<Numeric>(i + m + 1), static_cast<Numeric>(-2 * m));
    }

    for (Index i = m; i < NLeg; i++) {
      for (Index j = 0; j < N; j++) {
        asso_leg_term_pos(i - m, j) = Legendre::assoc_legendre(i, m, mu_arr[j]);
        asso_leg_term_neg(i - m, j) =
            asso_leg_term_pos(i - m, j) * ((i - m) % 2 ? -1.0 : 1.0);
      }
      asso_leg_term_mu0[i - m] = Legendre::assoc_legendre(i, m, -mu0);
    }

    const bool all_asso_leg_term_pos_finite =
        std::all_of(asso_leg_term_pos.elem_begin(),
                    asso_leg_term_pos.elem_end(),
                    [](auto& x) { return std::isfinite(x); });

    TIMEMACRO(dtstart += TimeStep(Time{} - tm).count());

    for (Index l = 0; l < NLayers; l++) {
      TIMEMACRO(Time lstart{});
      for (Index i = 0; i < NLeg - m; i++) {
        weighted_asso_Leg_coeffs_l[i] =
            weighted_scaled_Leg_coeffs(l, i + m) * fac[i];
      }

      const Numeric scaled_omega_l = scaled_omega_arr[l];

      TIMEMACRO(dtlstart += TimeStep(Time{} - lstart).count());

      if (all_asso_leg_term_pos_finite and
          std::any_of(weighted_asso_Leg_coeffs_l.elem_begin(),
                      weighted_asso_Leg_coeffs_l.elem_end(),
                      Cmp::gt(0))) {
        TIMEMACRO(Time D{});
        einsum<"ij", "j", "ji">(
            D_temp, weighted_asso_Leg_coeffs_l, asso_leg_term_pos);
        mult(D_pos, D_temp, asso_leg_term_pos, 0.5 * scaled_omega_l);
        mult(D_neg, D_temp, asso_leg_term_neg, 0.5 * scaled_omega_l);
        TIMEMACRO(dtD += TimeStep(Time{} - D).count());

        TIMEMACRO(Time rrsqr{});
        einsum<"ij", "i", "ij", "j">(sqr, inv_mu_arr.slice(0, N), D_neg, W);
        einsum<"ij", "i", "ij", "j">(apb, inv_mu_arr.slice(0, N), D_pos, W);
        apb.diagonal() -= inv_mu_arr.slice(0, N);

        auto K = Km[l];
        auto G = Gm[l];

        amb = apb;   // still just alpha
        apb += sqr;  // sqr is beta
        amb -= sqr;
        mult(sqr, amb, apb);
        TIMEMACRO(dtsqr += TimeStep(Time{} - rrsqr).count());

        //FIXME: The matrix produces real eigen values, a specialized solver might be good
        TIMEMACRO(Time diag{});
        ::diagonalize_inplace(
            amb, K.slice(0, N), K.slice(N, N), sqr, diag_work);
        TIMEMACRO(dtdiag += TimeStep(Time{} - diag).count());

        TIMEMACRO(Time collectG{});
        G(Range(0, N), Range(0, N)) = amb;
        for (Index i = 0; i < N; i++) {
          G[i].slice(0, N) *= 0.5;
          G[i].slice(N, N) = G[i].slice(0, N);
        }
        for (Index j = 0; j < N; j++) {
          const Numeric sqrt_x = std::sqrt(K[j]);
          K[j] = -sqrt_x;
          K[j + N] = sqrt_x;
        }

        mult(GmG, apb, G.slice(0, N));
        for (Index j = 0; j < NQuad; j++) {
          GmG(joker, j) /= K[j];
        }
        G.slice(N, N) = G.slice(0, N);
        G.slice(0, N) -= GmG;
        G.slice(N, N) += GmG;
        TIMEMACRO(dtcollectG += TimeStep(Time{} - collectG).count());

        TIMEMACRO(Time beam{});
        if (has_beam_source) {
          einsum<"i", "i", "i", "">(
              X_temp,
              weighted_asso_Leg_coeffs_l,
              asso_leg_term_mu0,
              (scaled_omega_l * I0 * (2 - m_equals_0_bool) /
               (4 * Constant::pi)));

          mult(jvec.slice(0, N).reshape_as(1, N), xtemp, asso_leg_term_pos, -1);
          jvec.slice(0, N) *= inv_mu_arr.slice(0, N);

          mult(jvec.slice(N, N).reshape_as(1, N), xtemp, asso_leg_term_neg);
          jvec.slice(N, N) *= inv_mu_arr.slice(0, N);

          std::copy(G.elem_begin(), G.elem_end(), Gml.elem_begin());
          solve_inplace(jvec, Gml, solve_work);

          for (Index j = 0; j < NQuad; j++) {
            jvec[j] /= (1.0 / mu0 + K[j]);
          }

          mult(Bm[l], G, jvec, -1);
        }
        TIMEMACRO(dtbeam += TimeStep(Time{} - beam).count());
      } else {
        TIMEMACRO(Time tspecial{});
        auto G = Gm[l];
        for (Index i = 0; i < N; i++) {
          G(i + N, i) = 1;
          G(i, i + N) = 1;
        }
        for (Index i = 0; i < NQuad; i++) {
          Km(l, i) = -1 / mu_arr[i];
        }
        TIMEMACRO(dtspecial += TimeStep(Time{} - tspecial).count());
      }
    }

    TIMEMACRO(dtm += TimeStep(Time{} - tm).count());
  }

  TIMEMACRO(const auto ratio = [](auto a, auto b) {
    return std::round(1000 * a / b) / 10;
  });
  TIMEMACRO(std::cout << "Total Diagonalize time: " << dtm << '\n');
  TIMEMACRO(std::cout << "Ratio in START:   " << ratio(dtstart, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in SPECIAL: " << ratio(dtspecial, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in LSTART:  " << ratio(dtlstart, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in BEAM:    " << ratio(dtbeam, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in D:       " << ratio(dtD, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in SQR:     " << ratio(dtsqr, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in DIAG:    " << ratio(dtdiag, dtm) << '%'
                      << '\n');
  TIMEMACRO(std::cout << "Ratio in COLLECT: " << ratio(dtcollectG, dtm) << '%'
                      << '\n');
}

/** Computes the IMS factors
 * 
 * Dependent on:
 * - omega_arr
 * - tau_arr
 * - f_arr
 * - mu0
 * - Leg_coeffs_all
 *
 * Modifies:
 * - scaled_mu0
 * - Leg_coeffs_residue_avg
 * - omega_avg
 * - f_avg
 */
void main_data::set_ims_factors() {
  const Numeric sum1 = omega_arr * tau_arr.vec();
  omega_avg = sum1 / sum(tau_arr);
  const Numeric sum2 = f_arr.empty() ? 0.0
                                     : einsum<Numeric, "", "i", "i", "i">(
                                           {}, f_arr, omega_arr, tau_arr);
  f_avg = sum2 / sum1;
  for (Index i = 0; i < NLeg_all; i++) {
    Numeric sum3 = 0.0;
    if (not f_arr.empty()) {
      if (i < NLeg) {
        for (Index j = 0; j < NLayers; j++) {
          sum3 += f_arr[j] * omega_arr[j] * tau_arr[j];
        }
      } else {
        for (Index j = 0; j < NLayers; j++) {
          sum3 += Leg_coeffs_all(j, i) * omega_arr[j] * tau_arr[j];
        }
      }
    }

    Leg_coeffs_residue_avg[i] = static_cast<Numeric>(2 * i + 1) *
                                (2 * sum3 / sum2 - Math::pow2(sum3 / sum2));
  }
  scaled_mu0 = mu0 / (1 - omega_avg * f_avg);
}

void main_data::set_scales() {
  std::transform(omega_arr.begin(),
                 omega_arr.end(),
                 f_arr.begin(),
                 scaled_omega_arr.begin(),
                 [](auto&& omega, auto&& f) { return 1.0 - omega * f; });

  for (Index i = 0; i < NLayers; i++) {
    scale_tau[i] = 1.0 - omega_arr[i] * f_arr[i];
  }

  scaled_tau_arr_with_0[0] = 0;
  einsum<"i", "i", "i">(
      scaled_tau_arr_with_0.slice(1, NLayers), tau_arr, scale_tau);

  for (Index i = 0; i < NLayers; i++) {
    for (Index j = 0; j < NLeg; j++) {
      weighted_scaled_Leg_coeffs(i, j) = static_cast<Numeric>(2 * j + 1) *
                                         (Leg_coeffs_all(i, j) - f_arr[i]) /
                                         (1 - f_arr[i]);
    }
  }

  for (Index i = 0; i < NLayers; i++) {
    scaled_omega_arr[i] = omega_arr[i] * (1.0 - f_arr[i]) / scale_tau[i];
  }
}

void main_data::set_weighted_Leg_coeffs_all() {
  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < NLeg_all; i++) {
      weighted_Leg_coeffs_all(joker, i) =
          static_cast<Numeric>(2 * i + 1) * Leg_coeffs_all(j, i);
    }
  }
}

void main_data::set_beam_source(const Numeric I0_) {
  has_beam_source = I0_ > 0;

  if (std::all_of(b_pos.elem_begin(), b_pos.elem_end(), Cmp::eq(0)) and
      not has_source_poly and has_beam_source) {
    I0_orig = I0_;
    I0 = 1;
  } else {
    I0_orig = 1;
    I0 = I0_;
  }
}

void main_data::check_input_size() const {
  ARTS_USER_ERROR_IF(
      tau_arr.size() != NLayers, tau_arr.size(), " vs ", NLayers);

  ARTS_USER_ERROR_IF(
      omega_arr.size() != NLayers, omega_arr.size(), " vs ", NLayers);

  ARTS_USER_ERROR_IF(
      (source_poly_coeffs.shape() != std::array{NLayers, Nscoeffs}),
      source_poly_coeffs.shape(),
      " vs (",
      NLayers,
      ", ",
      Nscoeffs,
      ')');

  ARTS_USER_ERROR_IF(
      f_arr.size() != NLayers, f_arr.size(), " vs ", NLayers, " or 0");

  ARTS_USER_ERROR_IF((Leg_coeffs_all.shape() != std::array{NLayers, NLeg_all}),
                     " vs (",
                     NLayers,
                     ", ",
                     NLeg_all,
                     ")");

  ARTS_USER_ERROR_IF((b_pos.shape() != std::array{NFourier, N}),
                     b_pos.shape(),
                     " vs (",
                     NFourier,
                     ", ",
                     N,
                     ')')

  ARTS_USER_ERROR_IF((b_neg.shape() != std::array{NFourier, N}),
                     b_neg.shape(),
                     " vs (",
                     NFourier,
                     ", ",
                     N,
                     ')')

  ARTS_USER_ERROR_IF(
      brdf_fourier_modes.size() != static_cast<std::size_t>(NBDRF),
      brdf_fourier_modes.size(),
      " vs ",
      NBDRF);
}

void main_data::check_input_value() const {
  ARTS_USER_ERROR_IF(tau_arr.front() <= 0.0,
                     "tau_arr must be strictly positive");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(omega_arr,
                          [](auto&& omega) { return omega >= 1 or omega < 0; }),
      "omega_arr must be [0, 1)");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Leg_coeffs_all,
                          [](auto&& x) {
                            return x[0] != 1 or
                                   std::ranges::any_of(x, [](auto&& u) {
                                     return std::abs<Numeric>(u) > 1;
                                   });
                          }),
      "Leg_coeffs_all must have 1 in the first column and be [-1, 1] elsewhere");

  ARTS_USER_ERROR_IF(I0 < 0, "I0 must be non-negative");

  ARTS_USER_ERROR_IF(phi0 < 0 or phi0 >= Constant::two_pi,
                     "phi0 must be [0, 2*pi)");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(f_arr, [](auto&& x) { return x > 1 or x < 0; }),
      "f_arr must be [0, 1]");

  ARTS_USER_ERROR_IF(mu0 < 0 or mu0 > 1, "mu0 must be [0, 1]");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          mu_arr.slice(0, N),
          [mu = mu0](auto&& x) { return std::abs(x - mu) < 1e-8; }),
      "mu0 in mu_arr, this creates a singularity.  Change NQuad or mu0.");
}

void main_data::update_all(const Numeric I0_) {
  check_input_value();

  set_weighted_Leg_coeffs_all();
  if (I0_ >= 0) set_beam_source(I0_);
  set_scales();
  set_ims_factors();
  diagonalize();
  solve_for_coefs();
}

main_data::main_data(const Index NLayers_,
                     const Index NQuad_,
                     const Index NLeg_,
                     const Index NFourier_,
                     const Index Nscoeffs_,
                     const Index NLeg_all_,
                     const Index NBDRF_)
    : NLayers(NLayers_),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(Nscoeffs_),
      NLeg_all(NLeg_all_),
      NBDRF(NBDRF_),
      has_source_poly(Nscoeffs > 0),
      is_multilayer(NLayers > 1),
      // User data
      tau_arr(NLayers),
      omega_arr(NLayers),
      f_arr(NLayers),
      source_poly_coeffs(NLayers, Nscoeffs),
      Leg_coeffs_all(NLayers, NLeg_all),
      b_pos(NFourier, N),
      b_neg(NFourier, N),
      brdf_fourier_modes(NBDRF),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLeg_all),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  Legendre::PositiveDoubleGaussLegendre(mu_arr.slice(0, N), W);

  std::transform(
      mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](auto&& x) {
        return -x;
      });
  std::transform(
      mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](auto&& x) {
        return 1.0 / x;
      });
}

main_data::main_data(const Index NQuad_,
                     const Index NLeg_,
                     const Index NFourier_,
                     AscendingGrid tau_arr_,
                     Vector omega_arr_,
                     Matrix Leg_coeffs_all_,
                     Matrix b_pos_,
                     Matrix b_neg_,
                     Vector f_arr_,
                     Matrix source_poly_coeffs_,
                     std::vector<BDRF> brdf_fourier_modes_,
                     Numeric mu0_,
                     Numeric I0_,
                     Numeric phi0_)
    : NLayers(tau_arr_.size()),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(source_poly_coeffs_.ncols()),
      NLeg_all(Leg_coeffs_all_.ncols()),
      NBDRF(brdf_fourier_modes_.size()),
      has_source_poly(Nscoeffs > 0),
      is_multilayer(NLayers > 1),
      has_beam_source(I0_ > 0),
      // User data
      tau_arr(std::move(tau_arr_)),
      omega_arr(std::move(omega_arr_)),
      f_arr(std::move(f_arr_)),
      source_poly_coeffs(std::move(source_poly_coeffs_)),
      Leg_coeffs_all(std::move(Leg_coeffs_all_)),
      b_pos(std::move(b_pos_)),
      b_neg(std::move(b_neg_)),
      brdf_fourier_modes(std::move(brdf_fourier_modes_)),
      mu0(mu0_),
      I0(I0_),
      phi0(phi0_),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLeg_all),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  Legendre::PositiveDoubleGaussLegendre(mu_arr.slice(0, N), W);

  std::transform(
      mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](auto&& x) {
        return -x;
      });
  std::transform(
      mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](auto&& x) {
        return 1.0 / x;
      });

  check_input_size();
  update_all(I0_);
}

[[nodiscard]] Index main_data::tau_index(const Numeric tau) const {
  const Index l =
      std::distance(tau_arr.begin(), std::ranges::lower_bound(tau_arr, tau));
  ARTS_USER_ERROR_IF(
      l == NLayers, "tau (", tau, ") must be at most ", tau_arr.back());
  return l;
}

void main_data::u(u_data& data, const Numeric tau, const Numeric phi) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NFourier, NQuad);
  for (Index i = 0; i < NFourier; i++) {
    for (Index j = 0; j < N; j++) {
      data.exponent(i, j) =
          std::exp(K_collect(i, l, j) * (scaled_tau - scaled_tau_arr_lm1));
      data.exponent(i, j + N) =
          std::exp(K_collect(i, l, j + N) * (scaled_tau - scaled_tau_arr_l));
    }
  }

  data.um.resize(NFourier, NQuad);
  static_assert(
      matpack::einsum_optpath<"mi", "mij", "mj">(),
      "On Failure, the einsum has been changed to not use optimal path");
  einsum<"mi", "mij", "mj">(
      data.um, GC_collect(joker, l, joker, joker), data.exponent);

  if (has_beam_source) {
    for (Index m = 0; m < NFourier; m++) {
      for (Index i = 0; i < NQuad; i++) {
        data.um(m, i) += std::exp(-scaled_tau / mu0) * B_collect(m, l, i);
      }
    }
  }

  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.um[0],
              data.src,
              tau,
              source_poly_coeffs[l],
              G_collect[0][l],
              K_collect[0][l],
              inv_mu_arr,
              0,
              1.0,
              1.0);
  }

  data.intensities.resize(NQuad);
  data.intensities = 0.0;
  for (Index m = 0; m < NFourier; m++) {
    const Numeric cp = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const auto umm = data.um[m];
    for (Index i = 0; i < NQuad; i++) {
      data.intensities[i] += umm[i] * cp;
    }
  }

  data.intensities *= I0_orig;
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NQuad);
  for (Index j = 0; j < N; j++) {
    data.exponent[j] =
        std::exp(K_collect(0, l, j) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[j + N] =
        std::exp(K_collect(0, l, j + N) * (scaled_tau - scaled_tau_arr_l));
  }

  data.u0.resize(NQuad);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0,
              data.src,
              tau,
              source_poly_coeffs[l],
              G_collect[0][l],
              K_collect[0][l],
              inv_mu_arr);
  } else {
    data.u0 = 0.0;
  }

  mult(data.u0, GC_collect(0, l, joker, joker), data.exponent, 1.0, 1.0);

  if (has_beam_source) {
    const auto tmp = B_collect(0, l, joker);
    std::transform(tmp.elem_begin(),
                   tmp.elem_end(),
                   data.u0.elem_begin(),
                   data.u0.elem_begin(),
                   [scl = std::exp(-scaled_tau / mu0)](auto&& x, auto&& y) {
                     return scl * x + y;
                   });
  }

  data.u0 *= I0_orig;
}

Numeric calculate_nu(const Numeric mu,
                     const Numeric phi,
                     const Numeric mu_p,
                     const Numeric phi_p) {
  const Numeric scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi);
  return mu * mu_p + scl * std::sqrt(1.0 - mu * mu);
}

void calculate_nu(Vector& nu,
                  const ExhaustiveConstVectorView& mu,
                  const Numeric phi,
                  const Numeric mu_p,
                  const Numeric phi_p) {
  nu.resize(mu.size());

  std::transform(
      mu.begin(),
      mu.end(),
      nu.begin(),
      [mu_p, scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi)](
          auto&& x) { return x * mu_p + scl * std::sqrt(1.0 - x * x); });
}

void main_data::TMS(tms_data& data,
                    const Numeric tau,
                    const Numeric phi) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  // mathscr_B
  calculate_nu(data.nu, mu_arr, phi, -mu0, phi0);

  data.mathscr_B.resize(NLayers, NQuad);
  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < NQuad; i++) {
      const Numeric p_true =
          Legendre::legendre_sum(weighted_Leg_coeffs_all[j], data.nu[i]);
      const Numeric p_trun =
          Legendre::legendre_sum(weighted_scaled_Leg_coeffs[j], data.nu[i]);
      data.mathscr_B(j, i) = (scaled_omega_arr[j] * I0) / (4 * Constant::pi) *
                             (mu0 / (mu0 + mu_arr[i])) *
                             (p_true / (1.0 - f_arr[j]) - p_trun);
    }
  }

  data.TMS.resize(NQuad);
  const Numeric exptau = std::exp(-scaled_tau / mu0);
  for (Index i = 0; i < N; i++) {
    data.TMS[i] =
        exptau - std::exp((scaled_tau - scaled_tau_arr_l) / mu_arr[i] -
                          scaled_tau_arr_l / mu0);
    data.TMS[i + N] =
        exptau - std::exp((scaled_tau_arr_lm1 - scaled_tau) / mu_arr[i] -
                          scaled_tau_arr_lm1 / mu0);
  }
  data.TMS *= data.mathscr_B[l];

  if (tau_arr.size() > 1) {
    data.contribution_from_other_layers_pos.resize(N, NLayers);
    data.contribution_from_other_layers_pos = 0;
    data.contribution_from_other_layers_neg.resize(N, NLayers);
    data.contribution_from_other_layers_neg = 0;
    for (Index i = 0; i < NLayers; i++) {
      if (l > i) {
        // neg
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_neg(j, i) =
              (data.mathscr_B(i, j + N) *
               (std::exp((scaled_tau_arr_with_0[i] - scaled_tau) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau_arr_with_0[i] - scaled_tau) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0)));
        }
      } else if (l < i) {
        // pos
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_pos(j, i) =
              (data.mathscr_B(i, j) *
               (std::exp((scaled_tau - scaled_tau_arr_with_0[i]) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau - scaled_tau_arr_with_0[i]) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0)));
        }
      } else {
        continue;
      }
    }

    for (Index i = 0; i < N; i++) {
      data.TMS[i + 0] += sum(data.contribution_from_other_layers_pos[i]);
      data.TMS[i + N] += sum(data.contribution_from_other_layers_neg[i]);
    }
  }
}

void main_data::IMS(Vector& ims, const Numeric tau, const Numeric phi) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  ims.resize(N);
  for (Index i = 0; i < N; i++) {
    const Numeric nu = calculate_nu(mu_arr[i + N], phi, -mu0, phi0);
    const Numeric x = 1.0 / mu_arr[i] - 1.0 / scaled_mu0;
    const Numeric chi = (1 / (mu_arr[i] * scaled_mu0 * x)) *
                        ((tau - 1 / x) * std::exp(-tau / scaled_mu0) +
                         std::exp(-tau / mu_arr[i]) / x);
    ims[i] = (I0 / (4 * Constant::pi) * Math::pow2(omega_avg * f_avg) /
              (1 - omega_avg * f_avg) *
              Legendre::legendre_sum(Leg_coeffs_residue_avg, nu)) *
             chi;
  }
}

void main_data::u_corr(u_data& u_data,
                       Vector& ims,
                       tms_data& tms_data,
                       const Numeric tau,
                       const Numeric phi) const {
  TMS(tms_data, tau, phi);
  IMS(ims, tau, phi);
  u(u_data, tau, phi);

  for (Index i = 0; i < N; i++) {
    u_data.intensities[i] += I0_orig * tms_data.TMS[i];
  }

  for (Index i = N; i < NQuad; i++) {
    u_data.intensities[i] += I0_orig * (tms_data.TMS[i] + ims[i - N]);
  }
}

Numeric main_data::flux_up(flux_data& data, const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau > tau_arr.back(),
                     "tau (",
                     tau,
                     ") must be less than the last layer (",
                     tau_arr.back(),
                     ")");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.u0_pos.resize(N);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0_pos,
              data.src,
              tau,
              source_poly_coeffs[l],
              G_collect[0][l],
              K_collect[0][l],
              inv_mu_arr);
  } else {
    data.u0_pos = 0.0;
  }

  if (has_beam_source) {
    for (Index i = 0; i < N; i++) {
      data.u0_pos[i] += B_collect(0, l, i) * std::exp(-scaled_tau / mu0);
    }
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect(0, l, i) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect(0, l, i + N) * (scaled_tau - scaled_tau_arr_l));
  }

  mult(data.u0_pos,
       GC_collect(0, l, Range(0, N), joker),
       data.exponent,
       1.0,
       1.0);

  return Constant::two_pi * I0_orig *
         einsum<Numeric, "", "i", "i", "i">(
             {}, mu_arr.slice(0, N), W, data.u0_pos);
}

std::pair<Numeric, Numeric> main_data::flux_down(flux_data& data,
                                                 const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau > tau_arr.back(),
                     "tau (",
                     tau,
                     ") must be less than the last layer (",
                     tau_arr.back(),
                     ")");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.u0_neg.resize(N);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0_neg,
              data.src,
              tau,
              source_poly_coeffs[l],
              G_collect[0][l],
              K_collect[0][l],
              inv_mu_arr,
              N);
  } else {
    data.u0_neg = 0.0;
  }

  const Numeric direct_beam =
      has_beam_source ? I0 * mu0 * std::exp(-tau / mu0) : 0;
  const Numeric direct_beam_scaled =
      has_beam_source ? I0 * mu0 * std::exp(-scaled_tau / mu0) : 0;
  if (has_beam_source) {
    for (Index i = 0; i < N; i++) {
      data.u0_neg[i] += B_collect(0, l, i + N) * std::exp(-scaled_tau / mu0);
    }
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect(0, l, i) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect(0, l, i + N) * (scaled_tau - scaled_tau_arr_l));
  }

  mult(data.u0_neg,
       GC_collect(0, l, Range(N, N), joker),
       data.exponent,
       1.0,
       1.0);

  return {I0_orig *
              (Constant::two_pi * einsum<Numeric, "", "i", "i", "i">(
                                      {}, mu_arr.slice(0, N), W, data.u0_neg) -
               direct_beam + direct_beam_scaled),
          I0_orig * I0 * direct_beam};
}

void main_data::gridded_flux(ExhaustiveVectorView flux_up,
                             ExhaustiveVectorView flux_do,
                             ExhaustiveVectorView flux_dd) const {
  Vector u0(NQuad);
  Vector exponent(NQuad, 1);
  mathscr_v_data src(NQuad, Nscoeffs);

  for (Index l = 0; l < NLayers; l++) {
    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];

    if (has_source_poly) {
      mathscr_v(u0,
                src,
                tau_arr[l],
                source_poly_coeffs[l],
                G_collect[0][l],
                K_collect[0][l],
                inv_mu_arr);
    } else {
      u0 = 0.0;
    }

    const Numeric direct_beam =
        has_beam_source ? I0 * mu0 * std::exp(-tau_arr[l] / mu0) : 0;
    const Numeric direct_beam_scaled =
        has_beam_source ? I0 * mu0 * std::exp(-scaled_tau_arr_l / mu0) : 0;
    if (has_beam_source) {
      for (Index i = 0; i < NQuad; i++) {
        u0[i] += B_collect(0, l, i) * std::exp(-scaled_tau_arr_l / mu0);
      }
    }

    for (Index i = 0; i < N; i++) {
      exponent[i] = std::exp(K_collect(0, l, i) *
                             (scaled_tau_arr_l - scaled_tau_arr_lm1));
    }

    mult(u0, GC_collect(0, l, joker, joker), exponent, 1.0, 1.0);

    flux_up[l] = Constant::two_pi * I0_orig *
                 einsum<Numeric, "", "i", "i", "i">(
                     {}, mu_arr.slice(0, N), W, u0.slice(0, N));
    flux_do[l] = I0_orig * (Constant::two_pi *
                                einsum<Numeric, "", "i", "i", "i">(
                                    {}, mu_arr.slice(0, N), W, u0.slice(N, N)) -
                            direct_beam + direct_beam_scaled);
    flux_dd[l] = I0_orig * I0 * direct_beam;
  }
}

void main_data::gridded_u(ExhaustiveTensor3View out, const Vector& phi) const {
  Matrix exponent(NFourier, NQuad, 1);
  Matrix um(NFourier, NQuad);
  mathscr_v_data src(NQuad, Nscoeffs);

  const Index Nphi = phi.size();
  Matrix cp(Nphi, NFourier);
  for (Index p = 0; p < phi.size(); p++) {
    for (Index m = 0; m < NFourier; m++) {
      cp(p, m) = I0_orig * std::cos(static_cast<Numeric>(m) * (phi0 - phi[p]));
    }
  }

  for (Index l = 0; l < NLayers; l++) {
    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];

    for (Index i = 0; i < NFourier; i++) {
      for (Index j = 0; j < N; j++) {
        exponent(i, j) = std::exp(K_collect(i, l, j) *
                                  (scaled_tau_arr_l - scaled_tau_arr_lm1));
      }
    }

    static_assert(
        matpack::einsum_optpath<"mi", "mij", "mj">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"mi", "mij", "mj">(um, GC_collect(joker, l, joker, joker), exponent);

    if (has_beam_source) {
      for (Index m = 0; m < NFourier; m++) {
        for (Index i = 0; i < NQuad; i++) {
          um(m, i) += std::exp(-scaled_tau_arr_l / mu0) * B_collect(m, l, i);
        }
      }
    }

    if (has_source_poly) {
      mathscr_v(um[0],
                src,
                tau_arr[l],
                source_poly_coeffs[l],
                G_collect[0][l],
                K_collect[0][l],
                inv_mu_arr,
                0,
                1.0,
                1.0);
    }

    static_assert(
        matpack::einsum_optpath<"pi", "im", "pm">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(out[l], transpose(um), cp);
  }
}

void main_data::ungridded_flux(ExhaustiveVectorView flux_up,
                               ExhaustiveVectorView flux_do,
                               ExhaustiveVectorView flux_dd,
                               const AscendingGrid& tau) const {
  ARTS_USER_ERROR_IF(tau.front() < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau.back() > tau_arr.back(),
                     "tau (",
                     tau.back(),
                     ") must be no more than the last layer (",
                     tau_arr.back(),
                     ")");

  Vector u0(NQuad);
  Vector exponent(NQuad, 1);
  mathscr_v_data src(NQuad, Nscoeffs);

  Index l = tau_index(tau.front());
  for (Index il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
    const Numeric scaled_tau =
        scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    if (has_source_poly) {
      mathscr_v(u0,
                src,
                tau[il],
                source_poly_coeffs[l],
                G_collect[0][l],
                K_collect[0][l],
                inv_mu_arr);
    } else {
      u0 = 0.0;
    }

    const Numeric direct_beam =
        has_beam_source ? I0 * mu0 * std::exp(-tau[il] / mu0) : 0;
    const Numeric direct_beam_scaled =
        has_beam_source ? I0 * mu0 * std::exp(-scaled_tau / mu0) : 0;
    if (has_beam_source) {
      for (Index i = 0; i < NQuad; i++) {
        u0[i] += B_collect(0, l, i) * std::exp(-scaled_tau / mu0);
      }
    }

    for (Index i = 0; i < N; i++) {
      exponent[i] =
          std::exp(K_collect(0, l, i) * (scaled_tau - scaled_tau_arr_lm1));
      exponent[i + N] =
          std::exp(K_collect(0, l, i + N) * (scaled_tau - scaled_tau_arr_l));
    }

    mult(u0, GC_collect(0, l, joker, joker), exponent, 1.0, 1.0);

    flux_up[il] = Constant::two_pi * I0_orig *
                  einsum<Numeric, "", "i", "i", "i">(
                      {}, mu_arr.slice(0, N), W, u0.slice(0, N));
    flux_do[il] =
        I0_orig *
        (Constant::two_pi * einsum<Numeric, "", "i", "i", "i">(
                                {}, mu_arr.slice(0, N), W, u0.slice(N, N)) -
         direct_beam + direct_beam_scaled);
    flux_dd[il] = I0_orig * I0 * direct_beam;
  }
}

void main_data::ungridded_u(ExhaustiveTensor3View out,
                            const AscendingGrid& tau,
                            const Vector& phi) const {
  ARTS_USER_ERROR_IF(tau.front() < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau.back() > tau_arr.back(),
                     "tau (",
                     tau.back(),
                     ") must be no more than the last layer (",
                     tau_arr.back(),
                     ")");

  Matrix exponent(NFourier, NQuad, 1);
  Matrix um(NFourier, NQuad);
  mathscr_v_data src(NQuad, Nscoeffs);

  const Index Nphi = phi.size();
  Matrix cp(Nphi, NFourier);
  for (Index p = 0; p < phi.size(); p++) {
    for (Index m = 0; m < NFourier; m++) {
      cp(p, m) = I0_orig * std::cos(static_cast<Numeric>(m) * (phi0 - phi[p]));
    }
  }

  Index l = tau_index(tau.front());
  for (Index il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
    const Numeric scaled_tau =
        scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    for (Index i = 0; i < NFourier; i++) {
      for (Index j = 0; j < N; j++) {
        exponent(i, j) =
            std::exp(K_collect(i, l, j) * (scaled_tau - scaled_tau_arr_lm1));
        exponent(i, j + N) =
            std::exp(K_collect(i, l, j + N) * (scaled_tau - scaled_tau_arr_l));
      }
    }

    static_assert(
        matpack::einsum_optpath<"mi", "mij", "mj">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"mi", "mij", "mj">(um, GC_collect(joker, l, joker, joker), exponent);

    if (has_beam_source) {
      for (Index m = 0; m < NFourier; m++) {
        for (Index i = 0; i < NQuad; i++) {
          um(m, i) += std::exp(-scaled_tau / mu0) * B_collect(m, l, i);
        }
      }
    }

    if (has_source_poly) {
      mathscr_v(um[0],
                src,
                tau[il],
                source_poly_coeffs[l],
                G_collect[0][l],
                K_collect[0][l],
                inv_mu_arr,
                0,
                1.0,
                1.0);
    }

    static_assert(
        matpack::einsum_optpath<"pi", "im", "pm">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(out[il], transpose(um), cp);
  }
}
}  // namespace disort