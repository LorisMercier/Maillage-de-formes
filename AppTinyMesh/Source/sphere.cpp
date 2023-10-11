// Sphere

// Self include
#include "sphere.h"

Sphere::Sphere(){
   c = Vector(0,0,0);
   r = 3;
}

/*!
 * \brief Create a sphere given a center point, the length of the radius.
 * \param _c Center.
 * \param _r radius lenght.
 */
Sphere::Sphere(const Vector& _c, double _r){
    c = _c;
    r = _r;
}


Sphere::~Sphere(){

}
