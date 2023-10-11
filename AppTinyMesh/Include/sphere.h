#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"


class Sphere
{
public:
    Vector c;//!< sphere center
    double r;//!< radius size
  //! Empty.
  Sphere();

  Sphere(const Vector&, double);

  //! Empty.
  ~Sphere();
};

