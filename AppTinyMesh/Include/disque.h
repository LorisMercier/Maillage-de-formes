#ifndef DISQUE_H
#define DISQUE_H

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"


class Disque
{
public:
    Vector c;//disk center
    double r;//radius size
  //! Empty.
  Disque();

  Disque(const Vector&, double);

  //! Empty.
  ~Disque();
};


#endif // DISQUE_H
