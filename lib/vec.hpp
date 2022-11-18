#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

using std::sqrt;

template <typename T, size_t Dim>
class vecTDim {
public:

  template <typename... Args>
  vecTDim(Args... args) {
    e = {static_cast<T>(args)...};
  }

  T x() const {
    return e[0];
  }
  T y() const {
    return e[1];
  }
  T z() const {
    return e[2];
  }
  T w() const {
    return e[3];
  }

  vecTDim operator-() const {
    vecTDim tmp = *this;
    return tmp *= -1;
  }
  T operator[](int i) const {
    return e[i];
  }
  T& operator[](int i) {
    return e[i];
  }

  vecTDim& operator+=(const vecTDim& v) {
    for (size_t i = 0; i < Dim; i++) {
      e[i] += v.e[i];
    }
    return *this;
  }

  vecTDim& operator*=(const T t) {
    for (auto& c : e) {
      c *= t;
    }
    return *this;
  }

  vecTDim& operator*=(const vecTDim& v) {
    for (size_t i = 0; i < Dim; i++) {
      e[i] *= v.e[i];
    }
    return *this;
  }

  vecTDim& operator/=(const T t) {
    return *this *= 1 / t;
  }

  T length() const {
    return sqrt(length_squared());
  }

  T length_squared() const {
    T sum = 0;
    for (auto c : e) {
      sum += c * c;
    }
    return sum;
  }

  /****************************************************************************
  *                              Static members                              *
  ****************************************************************************/

  static constexpr size_t DIM = Dim;

  static vecTDim unit() {
    vecTDim res;
    for (auto& x : res) {
      x = 1;
    }
    return res;
  }

  static vecTDim zero() {
    vecTDim res;
    for (auto& x : res) {
      x = 0;
    }
    return res;
  }

private:
  vecTDim() = default;

public:
  std::array<T, Dim> e;
};

/*******************************************************************************
 *                                Type aliases                                 *
 *******************************************************************************/

template<typename T>
using v2 = vecTDim<T, 2>;
template<typename T>
using v3 = vecTDim<T, 3>;

using v3d = v3<double>;
using v3f = v3<float>;
using v3i = v3<std::int32_t>;
using v3ui = v3<std::uint32_t>;
using v3l = v3<std::int64_t>;
using v3ul = v3<std::uint64_t>;

using v2d = v2<double>;
using v2f = v2<float>;
using v2i = v2<std::int32_t>;
using v2ui = v2<std::uint32_t>;
using v2l = v2<std::int64_t>;
using v2ul = v2<std::uint64_t>;

/*******************************************************************************
 *                              Utility functions                              *
 *******************************************************************************/

template <typename T, size_t Dim>
std::ostream& operator<<(std::ostream& out, const vecTDim<T, Dim>& v) {
  for (auto c : v.e) {
    out << c << ' ';
  }
  return out;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator+(const vecTDim<T, Dim>& u, const vecTDim<T, Dim>& v) {
  auto tmp(u);
  return tmp += v;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator-(const vecTDim<T, Dim>& u, const vecTDim<T, Dim>& v) {
  auto tmp(u);
  return tmp += -v;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator*(const vecTDim<T, Dim>& u, const vecTDim<T, Dim>& v) {
  auto tmp(u);
  return tmp *= v;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator*(T t, const vecTDim<T, Dim>& v) {
  auto tmp(v);
  return tmp *= t;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator*(const vecTDim<T, Dim>& v, double t) {
  return t * v;
}

template <typename T, size_t Dim>
vecTDim<T, Dim> operator/(vecTDim<T, Dim> v, double t) {
  return (1 / t) * v;
}

template <typename T, size_t Dim>
double dot(const vecTDim<T, Dim>& u, const vecTDim<T, Dim>& v) {
  T sum = 0;
  for (size_t i = 0; i < Dim; i++) {
    sum += u.e[i] * v.e[i];
  }
  return sum;
}

template <typename T>
vecTDim<T, 3> cross(const vecTDim<T, 3>& u, const vecTDim<T, 3>& v) {
  return vecTDim<T, 3>(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                       u.e[2] * v.e[0] - u.e[0] * v.e[2],
                       u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

template <typename T, size_t Dim>
vecTDim<T, Dim> unit_vector(vecTDim<T, Dim> v) {
  return v / v.length();
}

/*******************************************************************************
*                         Structured bindings support                         *
*******************************************************************************/

template <typename T, size_t Dim>
struct std::tuple_size<vecTDim<T, Dim>> {
  static constexpr size_t value = Dim;
};

template <size_t Idx, typename T, size_t Dim>
struct std::tuple_element<Idx, vecTDim<T, Dim>> {
  using type = T;
};

template <size_t Idx, typename T, size_t Dim>
T get(const vecTDim<T, Dim>& vec) {
  return vec.e[Idx];
}

template <size_t Idx, typename T, size_t Dim>
T& get(vecTDim<T, Dim>& vec) {
  return vec.e[Idx];
}
