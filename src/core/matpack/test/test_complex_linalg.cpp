#include <lin_alg.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace {
void require(const bool condition, const std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

void close(const Complex actual, const Complex expected, const std::string_view message) {
  require(std::abs(actual - expected) <= 2e-12 * std::max(1.0, std::abs(expected)), message);
}

void rejects(const std::function<void()>& operation, const std::string_view message) {
  try {
    operation();
  } catch (const std::exception& error) {
    require(not std::string_view(error.what()).empty(), "Validation error has no explanation");
    return;
  }
  throw std::runtime_error(std::string(message));
}

ComplexMatrix nonsymmetric() {
  ComplexMatrix A(3, 3);
  A[0, 0] = {3, 1};
  A[0, 1] = {2, -1};
  A[0, 2] = {-1, 2};
  A[1, 0] = {-1, 0};
  A[1, 1] = {4, -2};
  A[1, 2] = {0, 1};
  A[2, 0] = {2, 1};
  A[2, 1] = {-2, -1};
  A[2, 2] = {5, 3};
  return A;
}

ComplexVector product(StridedConstComplexMatrixView A, StridedConstComplexVectorView x) {
  ComplexVector out(A.nrows(), 0);
  for (Index i = 0; i < A.nrows(); ++i)
    for (Index j = 0; j < A.ncols(); ++j) out[i] += A[i, j] * x[j];
  return out;
}

void same(StridedConstComplexMatrixView A, StridedConstComplexMatrixView B) {
  require(A.shape() == B.shape(), "Changed matrix dimensions");
  for (Index i = 0; i < A.nrows(); ++i)
    for (Index j = 0; j < A.ncols(); ++j) require(A[i, j] == B[i, j], "Changed input matrix");
}

void eigen_residual(StridedConstComplexMatrixView A, StridedConstComplexMatrixView P, StridedConstComplexVectorView W) {
  for (Index j = 0; j < A.ncols(); ++j) {
    Numeric norm = 0;
    for (Index i = 0; i < A.nrows(); ++i) {
      Complex ap{};
      for (Index k = 0; k < A.ncols(); ++k) ap += A[i, k] * P[k, j];
      close(ap, P[i, j] * W[j], "Right-eigenvector residual is too large");
      norm += std::norm(P[i, j]);
    }
    close(norm, 1, "Eigenvector is not normalized");
  }
}

void ordinary() {
  const auto          A = nonsymmetric();
  const ComplexVector expected{{1, 2}, {-2, 1}, {3, -1}};
  const auto          b        = product(A, expected);
  const auto          original = A;
  ComplexVector       x(3);
  const Numeric       rcond = solve(x, A, b, 1e-12);
  require(rcond > 0 and rcond <= 1, "Invalid reciprocal condition estimate");
  for (Size i = 0; i < x.size(); ++i) close(x[i], expected[i], "Incorrect complex solution");
  same(A, original);

  auto rhs = b;
  solve(rhs, A, rhs);
  for (Size i = 0; i < rhs.size(); ++i) close(rhs[i], expected[i], "Aliased solve is incorrect");

  ComplexMatrix                P(3, 3);
  ComplexVector                W(3);
  complex_diagonalize_workdata work(3);
  diagonalize(P, W, A, work);
  eigen_residual(A, P, W);
  same(A, original);
  require(work.work.size() >= 6, "Eigenvalue workspace query was not retained");
  diagonalize(P, W, A, work);
  eigen_residual(A, P, W);

  auto alias = A;
  diagonalize(alias, W, alias);
  eigen_residual(A, alias, W);
}

void strided() {
  constexpr Complex   guard{17, -31};
  ComplexMatrix       backing(6, 7, guard), vectors(6, 7, guard);
  ComplexVector       right(7, guard), result(7, guard), values(7, guard);
  auto                A = backing[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto                P = vectors[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto                b = right[StridedRange(1, 3, 2)];
  auto                x = result[StridedRange(1, 3, 2)];
  auto                W = values[StridedRange(1, 3, 2)];
  const ComplexVector expected{{1, -2}, {2, 1}, {-1, 3}};
  A                   = nonsymmetric();
  b                   = product(A, expected);
  const auto original = backing;
  solve(x, A, b);
  for (Size i = 0; i < x.size(); ++i) close(x[i], expected[i], "Incorrect strided solution");
  diagonalize(P, W, A);
  eigen_residual(A, P, W);
  same(backing, original);
  for (Index i = 0; i < 7; i += 2) {
    require(result[i] == guard, "Strided solve overwrote an unselected value");
    require(values[i] == guard, "Strided eigenvalues overwrote an unselected value");
  }
  for (Index i = 0; i < 6; ++i)
    for (Index j = 0; j < 7; ++j)
      if (i % 2 != 0 or j % 2 == 0) require(vectors[i, j] == guard, "Strided eigenvectors overwrote padding");
}

void invalid_and_empty() {
  ComplexMatrix singular(2, 2, 1);
  ComplexVector b(2, 1), x(2, Complex{17, -31});
  rejects([&] { solve(x, singular, b); }, "Accepted a singular complex matrix");
  for (Complex value : x) require(value == Complex{17, -31}, "Failed solve changed its output");

  ComplexMatrix ill(2, 2, 0);
  ill[0, 0] = 1;
  ill[1, 1] = 1e-14;
  rejects([&] { solve(x, ill, b, 1e-12); }, "Accepted an ill-conditioned complex matrix");
  const Numeric rcond = solve(x, ill, b);
  require(std::abs(rcond / 1e-14 - 1) < 1e-12, "Incorrect condition estimate for diagonal matrix");
  close(x[1], 1e14, "Incorrect unthresholded solution");

  for (const Numeric threshold : {-1.0, 1.1, std::numeric_limits<Numeric>::quiet_NaN()})
    rejects([&] { solve(x, ill, b, threshold); }, "Accepted an invalid condition threshold");

  ComplexMatrix P(2, 2, 7), rectangular(2, 3);
  ComplexVector W(2, 7), short_vector(1);
  rejects([&] { solve(x, rectangular, b); }, "Accepted a nonsquare solve matrix");
  rejects([&] { solve(short_vector, ill, b); }, "Accepted incorrect solution dimensions");
  rejects([&] { diagonalize(P, W, rectangular); }, "Accepted a nonsquare eigenproblem");
  rejects([&] { diagonalize(P, short_vector, ill); }, "Accepted incorrect eigenvalue dimensions");
  complex_diagonalize_workdata wrong_work(1);
  rejects([&] { diagonalize(P, W, ill, wrong_work); }, "Accepted incorrect workspace dimensions");
  rejects([] { complex_diagonalize_workdata invalid(-1); }, "Accepted a negative workspace size");

  ill[0, 0] = {1, std::numeric_limits<Numeric>::infinity()};
  rejects([&] { solve(x, ill, b); }, "Accepted nonfinite solve matrix input");
  rejects([&] { diagonalize(P, W, ill); }, "Accepted nonfinite eigenproblem input");
  for (Complex value : W) require(value == Complex{7, 0}, "Failed diagonalization changed eigenvalues");
  for (Index i = 0; i < 2; ++i)
    for (Index j = 0; j < 2; ++j) require(P[i, j] == Complex{7, 0}, "Failed diagonalization changed eigenvectors");
  ill[0, 0] = 1;
  b[0]      = {std::numeric_limits<Numeric>::quiet_NaN(), 0};
  rejects([&] { solve(x, ill, b); }, "Accepted nonfinite right-hand side input");

  // A valid eigenproblem can have a singular eigenvector basis. The caller's
  // condition threshold must reject it before computing equivalent strengths.
  ComplexMatrix jordan(2, 2, 0);
  jordan[0, 0] = jordan[1, 1] = jordan[0, 1] = 1;
  diagonalize(P, W, jordan);
  b = 1;
  rejects([&] { solve(x, P, b, 1e-12); }, "Accepted a defective eigenvector basis");

  ComplexMatrix empty(0, 0);
  ComplexVector empty_vector(0);
  require(solve(empty_vector, empty, empty_vector) == 1, "Empty system has incorrect condition estimate");
  diagonalize(empty, empty_vector, empty);
  complex_diagonalize_workdata empty_work;
  diagonalize(empty, empty_vector, empty, empty_work);
}
}  // namespace

int main() try {
  ordinary();
  strided();
  invalid_and_empty();
  std::cout << "Complex LAPACK solve and diagonalization checks passed.\n";
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
