#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"


class Capsule
{
public:
  Vector c;//!< capsule center
  double r,h;//!< radius and height size
  //! Empty.
  Capsule();

  Capsule(const Vector&, double,double);

  //! Empty.
  ~Capsule();
};

