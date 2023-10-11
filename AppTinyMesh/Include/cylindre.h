#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"


class Cylindre
{
public:
    Vector c;//!< cylinder center
    double r,h;//!< radius and height size
  //! Empty.
  Cylindre();

  Cylindre(const Vector&, double,double);

  //! Empty.
  ~Cylindre();
};
