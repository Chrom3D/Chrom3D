#ifndef GUARD_COORDINATES_H
#define GUARD_COORDINATES_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <assert.h>

// Small vector/matrix types for bead coordinates. They replace the boost::numeric::ublas types
// used in earlier versions and provide only the operations Chrom3D needs (with the same names:
// v(i), m(i,j), size1() = number of rows, size2() = number of columns).

class CoordinateVector {
 public:
  explicit CoordinateVector(std::size_t n=0) : values(n, 0.0) {}
  double& operator()(std::size_t i) { return values[i]; }
  double operator()(std::size_t i) const { return values[i]; }
  std::size_t size() const { return values.size(); }
 private:
  std::vector<double> values;
};

class CoordinateMatrix {
 public:
  CoordinateMatrix(std::size_t rows=0, std::size_t cols=0) : nRows(rows), nCols(cols), values(rows*cols, 0.0) {}
  double& operator()(std::size_t i, std::size_t j) { return values[i*nCols + j]; }
  double operator()(std::size_t i, std::size_t j) const { return values[i*nCols + j]; }
  std::size_t size1() const { return nRows; }
  std::size_t size2() const { return nCols; }
  CoordinateVector row(std::size_t i) const {
    CoordinateVector res(nCols);
    for(std::size_t j=0; j!=nCols; j++) {
      res(j) = (*this)(i,j);
    }
    return res;
  }
 private:
  std::size_t nRows;
  std::size_t nCols;
  std::vector<double> values; // row-major
};

inline CoordinateVector operator-(const CoordinateVector& v1, const CoordinateVector& v2) {
  assert(v1.size() == v2.size());
  CoordinateVector res(v1.size());
  for(std::size_t i=0; i!=v1.size(); i++) {
    res(i) = v1(i) - v2(i);
  }
  return res;
}

inline CoordinateVector operator/(const CoordinateVector& v, double d) {
  CoordinateVector res(v.size());
  for(std::size_t i=0; i!=v.size(); i++) {
    res(i) = v(i) / d;
  }
  return res;
}

inline CoordinateMatrix operator+(const CoordinateMatrix& m1, const CoordinateMatrix& m2) {
  assert(m1.size1() == m2.size1() and m1.size2() == m2.size2());
  CoordinateMatrix res(m1.size1(), m1.size2());
  for(std::size_t i=0; i!=m1.size1(); i++) {
    for(std::size_t j=0; j!=m1.size2(); j++) {
      res(i,j) = m1(i,j) + m2(i,j);
    }
  }
  return res;
}

inline CoordinateMatrix operator-(const CoordinateMatrix& m1, const CoordinateMatrix& m2) {
  assert(m1.size1() == m2.size1() and m1.size2() == m2.size2());
  CoordinateMatrix res(m1.size1(), m1.size2());
  for(std::size_t i=0; i!=m1.size1(); i++) {
    for(std::size_t j=0; j!=m1.size2(); j++) {
      res(i,j) = m1(i,j) - m2(i,j);
    }
  }
  return res;
}

inline std::ostream& operator<<(std::ostream& os, const CoordinateMatrix& m) {
  for(std::size_t i=0; i!=m.size1(); i++) {
    for(std::size_t j=0; j!=m.size2(); j++) {
      os << (j == 0 ? "" : "\t") << m(i,j);
    }
    os << std::endl;
  }
  return os;
}

#endif
