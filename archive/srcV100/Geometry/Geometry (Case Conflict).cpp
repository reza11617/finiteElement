#include "geometry.h"

Geometry::Geometry(double x, double y, unsigned int m, unsigned int n) {
  dimentionX_ = x;
  dimentionY_ = y;
  NOEx_ = m;
  NOEy_ = n;
};

void Geometry::nodeDimentionVector(std::vector<double>& x, std::vector<double>& y) {
  // this function numbers the nodes the use method for numbering is
  /*
    6 --- 7 --- 8
    3 --- 4 --- 5
    0 --- 1 --- 2
  */
  double incr_x = dimentionX_/NOEx_;
  double incr_y = dimentionY_/NOEy_;
  for (unsigned int j = 0; j <= NOEy_; j++){
    for (unsigned int i = 0; i <= NOEx_; i++){
      x.push_back(i*incr_x);
      y.push_back(j*incr_y);
    }
  }
};

void Geometry::meshCantilever(std::vector<unsigned int>& elements)
{
  // this function meshes a contilever
  /*
    6 --- 7 --- 8
    |  3  |  4  |
    3 --- 4 --- 5
    |  1  |  2  |
    0 --- 1 --- 2
   */
  unsigned int x = NOEx_+1;
  unsigned int y = NOEy_;
  for (unsigned int i = 0; i < NOEy_; i++){
    for (unsigned int j = 0; j < NOEx_; j++){
      elements.push_back(i*x+j);
      elements.push_back(i*x+1+j);
      elements.push_back((i+1)*x+1+j);
      elements.push_back((i+1)*x+j);
    }
  }
};

