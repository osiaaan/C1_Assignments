#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>


class Vector : public std::vector<double> {
public:
  typedef std::vector<double> Base;

  Vector() : Base() {}
  Vector(int N) : Base(N) {}
  Vector(int N, double v) : Base(N, v) {}

  Vector operator+(const Vector& other) const;
  Vector operator-(const Vector& other) const;
  double operator*(const Vector& other) const;
  Vector operator*(double v) const;

  double norm() const;
  double maxNorm() const;


  void toFile(const std::string& filename, double L = 1.0, double tau = 0.0);
};

std::ostream& operator<<(std::ostream& s, const Vector& v);

#endif
