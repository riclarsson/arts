/**
 * @file wigner_functions.cc
 * @author Richard Larsson
 * @date 2013-06-19
 * 
 * @brief Wigner symbol interactions
 */

#include "wigner_functions.h"

#include <debug.h>
#include <wigxjpf_config.h>

#include <algorithm>

#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

void arts_wigner_thread_init(int max_two_j [[maybe_unused]]) {
#ifdef WIGXJPF_HAVE_THREAD
  wig_thread_temp_init(max_two_j);
#else
  wig_temp_init(max_two_j);
#endif
}

void arts_wigner_thread_free() { wig_temp_free(); }

int wigner_init_size(std::span<const Rational> x) noexcept {
  if (x.size() == 0) return 0;
  Rational out = x[0];
  for (Size i = 1; i < x.size(); i++) out = maxr(out, x[i]);
  return 1 + 2 * out.toInt();
}

Numeric wigner3j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational m1,
                 const Rational m2,
                 const Rational m3) {
  errno = 0;

  const int a = (2 * j1).toInt(), b = (2 * j2).toInt(), c = (2 * j3).toInt(),
            d = (2 * m1).toInt(), e = (2 * m2).toInt(), f = (2 * m3).toInt();
  double g;
  const int j = std::max({std::abs(a),
                          std::abs(b),
                          std::abs(c),
                          std::abs(d),
                          std::abs(e),
                          std::abs(f)}) *
                    3 / 2 +
                1;

  arts_wigner_thread_init(j);
  g = WIGNER3(a, b, c, d, e, f);
  arts_wigner_thread_free();

  if (errno == EDOM) {
    errno = 0;
    ARTS_USER_ERROR("Bad state, perhaps you need to call WignerInit?")
  }

  return Numeric(g);
}

Numeric wigner6j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational l1,
                 const Rational l2,
                 const Rational l3) {
  errno = 0;

  const int a = (2 * j1).toInt(), b = (2 * j2).toInt(), c = (2 * j3).toInt(),
            d = (2 * l1).toInt(), e = (2 * l2).toInt(), f = (2 * l3).toInt();
  double g;
  const int j = std::max({std::abs(a),
                          std::abs(b),
                          std::abs(c),
                          std::abs(d),
                          std::abs(e),
                          std::abs(f)});

  arts_wigner_thread_init(j);
  g = WIGNER6(a, b, c, d, e, f);
  arts_wigner_thread_free();

  if (errno == EDOM) {
    errno = 0;
    ARTS_USER_ERROR("Bad state, perhaps you need to call Wigner6Init?")
  }

  return Numeric(g);
}

std::pair<Rational, Rational> wigner_limits(std::pair<Rational, Rational> a,
                                            std::pair<Rational, Rational> b) {
  const bool invalid = a.first.isUndefined() or b.first.isUndefined();
  if (invalid) return {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};

  auto f = maxr(a.first, b.first);
  auto s = minr(a.second, b.second);

  if (f > s) return {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
  return {f, s};
}

Index make_wigner_ready(int largest, [[maybe_unused]] int fastest, int size) {
  if (size == 3) {
#if DO_FAST_WIGNER
    fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
#ifdef _OPENMP
    fastwigxj_thread_dyn_init(3, fastest);
#else
    fastwigxj_dyn_init(3, fastest);
#endif
#endif
    wig_table_init(largest, 3);

    return largest;
  }

  if (size == 6) {
#if DO_FAST_WIGNER
    fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
    fastwigxj_load(FAST_WIGNER_PATH_6J, 6, NULL);
#ifdef _OPENMP
    fastwigxj_thread_dyn_init(3, fastest);
    fastwigxj_thread_dyn_init(6, fastest);
#else
    fastwigxj_dyn_init(3, fastest);
    fastwigxj_dyn_init(6, fastest);
#endif
#endif
    wig_table_init(largest * 2, 6);

    return largest;
  }

  return 0;
}

namespace {
constexpr int wigner3_size(const Rational& J) { return 1 + 2 * J.toInt(6); }

constexpr int wigner6_size(const Rational& J) { return J.toInt(4) + 1; }

constexpr Rational wigner3_revere_size(const int j) {
  return Rational{(j - 1) / 2, 6};
}

constexpr Rational wigner6_revere_size(const int j) {
  return Rational{j - 1, 4};
}
}  // namespace

extern "C" int wigxjpf_max_prime_decomp;

bool is_wigner_ready(int j) { return not(j > wigxjpf_max_prime_decomp); }

bool is_wigner3_ready(const Rational& J) {
  const int test = wigner3_size(J);
  return is_wigner_ready(test);
}

bool is_wigner6_ready(const Rational& J) {
  const int test = wigner6_size(J);  // nb. J can be half-valued
  return is_wigner_ready(test);
}

namespace {
constexpr Index pow_negative_one(Index x) noexcept { return (x % 2) ? -1 : 1; }
}  // namespace

Numeric dwigner3j(Index M, Index J1, Index J2, Index J) {
  auto CJM = [](Index j, Index m) {
    Numeric cjm = 1.;
    for (Index I = 1; I <= j; I++) cjm *= (1. - .5 / static_cast<Numeric>(I));
    for (Index K = 1; K <= m; K++)
      cjm *= static_cast<Numeric>(j + 1 - K) / static_cast<Numeric>(j + K);
    return cjm;
  };

  Numeric GCM = 0.;
  if (J1 < 0 or J2 < 0 or J < 0) return GCM;
  Index JM = J1 + J2;
  Index J0 = std::abs(J1 - J2);
  if (J > JM or J < J0) return GCM;
  Index JS = std::max(J1, J2);
  Index JI = std::min(J1, J2);
  Index MA = std::abs(M);
  if (MA > JI) return GCM;
  Index UN    = 1 - 2 * (JS % 2);
  Index QM    = M + M;
  Numeric CG0 = 0.;
  GCM         = static_cast<Numeric>(UN) *
        std::sqrt(CJM(JI, MA) / CJM(JS, MA) * CJM(J0, 0) /
                  static_cast<Numeric>(JS + JS + 1));
  Index AJ0    = J0;
  Index AJM    = JM + 1;
  Index AJ02   = AJ0 * AJ0;
  Index AJM2   = AJM * AJM;
  Numeric ACG0 = 0.;
  for (Index I = J0 + 1; I <= J; I++) {
    Index AI    = I;
    Index AI2   = AI * AI;
    Numeric ACG = std::sqrt((AJM2 - AI2) * (AI2 - AJ02));
    Numeric CG1 =
        (static_cast<Numeric>(QM) * static_cast<Numeric>(I + I - 1) * GCM -
         ACG0 * CG0) /
        ACG;
    CG0  = GCM;
    GCM  = CG1;
    ACG0 = ACG;
  }
  return GCM;
}

Numeric dwigner6j(Index A, Index B, Index C, Index D, Index F) {
  Numeric SIXJ, TERM;
  if (std::abs(A - C) > F or std::abs(B - D) > F or (A + C) < F or (B + D < F))
    goto x1000;
  switch (C - D + 2) {
    case 2:  goto x2;
    case 3:  goto x3;
    case 4:  goto x1000;
    default: goto x1;
  }

x1:
  switch (A - B + 2) {
    case 2:  goto x11;
    case 3:  goto x12;
    default: goto x10;
  }

x10:
  TERM =
      (static_cast<Numeric>(F + B + D + 1) * static_cast<Numeric>(F + B + D) *
       static_cast<Numeric>(B + D - F) * static_cast<Numeric>(B + D - (1 + F)));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D - 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C + F)) * std::sqrt(TERM);
  return SIXJ;

x11:
  TERM  = (static_cast<Numeric>(F + B + D + 1) *
          static_cast<Numeric>(F + B - D + 1) *
          static_cast<Numeric>(F + D - B) * static_cast<Numeric>(B + D - F));
  TERM /= (4. * static_cast<Numeric>(B) * static_cast<Numeric>(2 * B + 1) *
           static_cast<Numeric>(B + 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(2 * D - 1));
  SIXJ =
      static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * std::sqrt(TERM);
  return SIXJ;

x12:
  TERM =
      (static_cast<Numeric>(F + D - B) * static_cast<Numeric>(F + D - B - 1) *
       static_cast<Numeric>(F + B - D + 2) *
       static_cast<Numeric>(F + B - D + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(2 * D - 1) *
           static_cast<Numeric>(D) * static_cast<Numeric>(2 * D + 1));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x2:
  switch (A - B + 2) {
    case 2:  goto x21;
    case 3:  goto x22;
    default: goto x20;
  }

x20:
  TERM =
      (static_cast<Numeric>(F + B + D + 1) * static_cast<Numeric>(D + B - F) *
       static_cast<Numeric>(F + B - D) * static_cast<Numeric>(F + D - B + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(D + 1));
  SIXJ =
      static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * std::sqrt(TERM);
  return SIXJ;

x21:
  TERM = (static_cast<Numeric>(B) * static_cast<Numeric>(B + 1) +
          static_cast<Numeric>(D) * static_cast<Numeric>(D + 1) -
          static_cast<Numeric>(F) * static_cast<Numeric>(F + 1));
  TERM /=
      std::sqrt(4. * static_cast<Numeric>(B) * static_cast<Numeric>(B + 1) *
                static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(D) *
                static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * TERM;
  return SIXJ;

x22:
  TERM =
      (static_cast<Numeric>(F + D + B + 2) *
       static_cast<Numeric>(F + B - D + 1) *
       static_cast<Numeric>(B + D - F + 1) * static_cast<Numeric>(F + D - B));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(D) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x3:
  switch (A - B + 2) {
    case 2:  goto x31;
    case 3:  goto x32;
    default: goto x30;
  }

x30:
  TERM =
      (static_cast<Numeric>(F + B - D) * static_cast<Numeric>(F + B - D - 1) *
       static_cast<Numeric>(F + D - B + 2) *
       static_cast<Numeric>(F + D - B + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D + 1) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(2 * D + 3));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x31:
  TERM =
      (static_cast<Numeric>(F + D + B + 2) *
       static_cast<Numeric>(B + D - F + 1) *
       static_cast<Numeric>(F + D - B + 1) * static_cast<Numeric>(F + B - D));
  TERM /= (4. * static_cast<Numeric>(B) * static_cast<Numeric>(2 * B + 1) *
           static_cast<Numeric>(B + 1) * static_cast<Numeric>(2 * D + 1) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 3));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x32:
  TERM  = (static_cast<Numeric>(F + D + B + 3) *
          static_cast<Numeric>(F + B + D + 2) *
          static_cast<Numeric>(B + D - F + 2) *
          static_cast<Numeric>(B + D - F + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(2 * D + 3) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ  = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x1000:
  SIXJ = 0.;
  return SIXJ;
}

int WignerInformation::largest = -1;
int WignerInformation::fastest = -1;
bool WignerInformation::threej = false;
bool WignerInformation::sixj   = false;
bool WignerInformation::init   = false;

std::ostream& operator<<(std::ostream& os, const WignerInformation& wi) {
  os << "Wigner Information: ";
  if (wi.init) {
    os << "Initialized.  Run *WignerUnload* to clear initialization.";
    os << "\n  Max allowed 3j symbol: J = ";
    if (wi.threej) {
      os << wigner3_revere_size(wi.largest);
    } else {
      os << "N/A";
    }
    os << "\n  Max allowed 6j symbol: J = ";
    if (wi.sixj) {
      os << wigner6_revere_size(wi.largest);
    } else {
      os << "N/A";
    }

#if DO_FAST_WIGNER
    if (wi.fastest) {
      os << "\n  Fast symbols are active (up to " << wi.fastest << " symbols)";
    } else
#endif
    {
      os << "\n  Fast symbols are inactive";
    }
  } else {
    os << "Not initialized.  Run *WignerInit* first with appropriate sizes";
  }
  return os;
}

void WignerInformation::initalize() {
  ARTS_USER_ERROR_IF(init, "Must not be initialized.")

  if (sixj) {
    make_wigner_ready(largest, fastest, 6);
  } else if (threej) {
    make_wigner_ready(largest, fastest, 3);
  }

  init = true;
}

void WignerInformation::unload() {
  if (sixj) {
#if DO_FAST_WIGNER
    fastwigxj_unload(3);
    fastwigxj_unload(6);
#endif
    wig_table_free();
  } else if (threej) {
#if DO_FAST_WIGNER
    fastwigxj_unload(3);
#endif
    wig_table_free();
  }

  init = false;
}

WignerInformation::WignerInformation(int largest_symbol,
                                     int fastest_symbol,
                                     bool three,
                                     bool six) {
  if (init) {
    if (largest_symbol != largest or fastest_symbol != fastest or
        three != threej or six != sixj) {
      unload();
    } else {
      return;
    }
  }

  ARTS_USER_ERROR_IF(largest_symbol < 0,
                     "You must specify a non-negative integer for largest.");
  ARTS_USER_ERROR_IF(fastest_symbol < 0,
                     "You must specify a non-negative integer for fastest.");
  ARTS_USER_ERROR_IF(not(three or six),
                     "You must specify at least one of threej or sixj.");

  largest = largest_symbol;
  fastest = fastest_symbol;
  threej  = three;
  sixj    = six;

  initalize();
}

void WignerInformation::assert_valid_wigner3(const Rational J) {
  ARTS_USER_ERROR_IF(not init, "Must not be initialized.")

  ARTS_USER_ERROR_IF(
      not(threej or is_wigner3_ready(J)),
      "Wigner library not ready for Wigner 3j symbols with J = {}"
      "\nPlease initialize it properly.",
      J);
}

void WignerInformation::assert_valid_wigner6(const Rational J) {
  ARTS_USER_ERROR_IF(not init, "Must not be initialized.")

  ARTS_USER_ERROR_IF(
      not(sixj or is_wigner6_ready(J)),
      "Wigner library not ready for Wigner 6j symbols with J = {}"
      "\nPlease initialize it properly.",
      J);
}
