#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"


class Tore
{
public:
    Vector c;//!< torus center
    double r;//!< size of the radius from the center of the torus to the center of the tube
    double rtube;//!< outer tube radius size
  //! Empty.
  Tore();

  Tore(const Vector&, double,double);

  //! Empty.
  ~Tore();
};
