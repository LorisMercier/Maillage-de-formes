// Tore

// Self include
#include "tore.h"

Tore::Tore(){
   c = Vector(0,0,0);
   r = 10;
   rtube = 10;
}

/*!
 * \brief Create a torus given a center point, the length of the radius.
 * \param _c Center.
 * \param _r size of the radius from the center of the torus to the center of the tube.
 * \param _r2 outer tube radius size.
 */
Tore::Tore(const Vector& _c, double _r,double _r2){
    c = _c;
    r = _r;
    rtube = _r2;
}


Tore::~Tore(){

}
