#ifndef RECORDER_H
#define RECORDER_H

#include <iostream>
#include <string>
#include <fstream>
#include "../Sparse/Sparse.h"

class Recorder {


public:
  const void matrix(const std::string& , float* , unsigned int) const;
  const void SparseMatrix(const std::string& fileName , const Sparse& sp) const;
  static Recorder& File();
  
};

#endif 
