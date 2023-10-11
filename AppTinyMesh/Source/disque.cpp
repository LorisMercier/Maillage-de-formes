// Disque

// Self include
#include "disque.h"

Disque::Disque(){
   c = Vector(0,0,0);
   r = 10;
}

/*!
 * \brief Create a disc given a center point, the length of the radius.
 * \param _c Center.
 * \param _r radius lenght.
 */
Disque::Disque(const Vector& _c, double _r){
    c = _c;
    r = _r;
}


Disque::~Disque(){

}
