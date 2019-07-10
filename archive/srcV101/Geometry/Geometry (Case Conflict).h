#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>

// Geometry class
class Geometry {
private:
  double dimentionX_ ;
  double dimentionY_ ;
  unsigned int NOEx_;
  unsigned int NOEy_;
public:
  Geometry(double, double, unsigned int, unsigned int);
  void nodeDimentionVector(std::vector<double>&, std::vector<double>&);
  void meshCantilever(std::vector<unsigned int>&);
};

#endif
