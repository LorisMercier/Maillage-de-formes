// Capsule

// Self include
#include "capsule.h"

Capsule::Capsule(){
   c = Vector(0,0,0);
   r = 2;
   h = 10;
}

/*!
 * \brief Create a capsule given a center point, the length of the radius and the height.
 * \param _c Center.
 * \param _r radius lenght.
 * \param _h height lenght.
 */
Capsule::Capsule(const Vector& _c, double _r,double _h){
    c = _c;
    r = _r;
    h = _h;
}


Capsule::~Capsule(){

}
