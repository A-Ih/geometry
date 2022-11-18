#pragma once

#include <array>
#include <cassert>
#include <iostream>
#include <utility>

template <typename T, size_t N, size_t M>
class Matrix {
public:
  Matrix(const Matrix&) = default;
  Matrix& operator=(const Matrix&) = default;

  Matrix(std::initializer_list<T> args) : e{args} {
    assert(e.size() == args.size());
  }

  T& operator()(size_t i, size_t j) {
    assert(i < N && j < M);
    return e[i * N + j];
  }

  const T& operator()(size_t i, size_t j) const {
    assert(i < N && j < M);
    return e[i * N + j];
  }

  Matrix operator-() const {
    Matrix tmp = *this;
    return tmp *= -1;
  }

  Matrix& operator+=(const Matrix& v) {
    for (size_t i = 0; i < e.size(); i++) {
      e[i] += v.e[i];
    }
    return *this;
  }

  Matrix& operator*=(const T t) {
    for (auto& c : e) {
      c *= t;
    }
    return *this;
  }

  Matrix& operator*=(const Matrix& v) {
    for (size_t i = 0; i < e.size(); i++) {
      e[i] *= v.e[i];
    }
    return *this;
  }

  Matrix& operator/=(const T t) {
    return *this *= 1 / t;
  }

  /****************************************************************************
  *                              Static members                              *
  ****************************************************************************/

  static constexpr size_t ROWS = N;
  static constexpr size_t COLS = M;

  static Matrix unit() {
    Matrix res = zero();
    for (size_t i = 0; i < std::min(N, M); i++) {
      res[i * N + i] = 1;
    }
    return res;
  }

  static Matrix zero() {
    Matrix res;
    for (auto& x : res) {
      x = 0;
    }
    return res;
  }

private:
  Matrix() = default;

  std::array<T, N * M> e;
};

/*******************************************************************************
 *                                Type aliases                                 *
 *******************************************************************************/

template<typename T>
using m3 = Matrix<T, 3, 3>;
template<typename T>
using m2 = Matrix<T, 2, 2>;

using m3d = m3<double>;
using m3f = m3<float>;
using m3i = m3<std::int32_t>;
using m3ui = m3<std::uint32_t>;
using m3l = m3<std::int64_t>;
using m3ul = m3<std::uint64_t>;

using m2d = m2<double>;
using m2f = m2<float>;
using m2i = m2<std::int32_t>;
using m2ui = m2<std::uint32_t>;
using m2l = m2<std::int64_t>;
using m2ul = m2<std::uint64_t>;

/*******************************************************************************
 *                              Utility functions                              *
 *******************************************************************************/

template <typename T, size_t N, size_t M>
std::ostream& operator<<(std::ostream& out, const Matrix<T, N, M>& v) {
  for (auto c : v.e) {
    out << c << ' ';
  }
  return out;
}

// Print a matrix in row-major order
template<typename T, size_t N, size_t M, size_t K>
Matrix<T, N, K> mul(const Matrix<T, N, M>& lhs, const Matrix<T, M, K>& rhs) {
  using MRes = Matrix<T, N, K>;
  MRes res = MRes::zero();
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < K; j++) {
      for (size_t k = 0; k < M; k++) {
        res(i, j) += lhs(i, k) * rhs(k, j);
      }
    }
  }
  return res;
}
