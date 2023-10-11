#include "heightField.h"

HeightField::HeightField(const QString& filename,const int& echelle):echelle(echelle){
    img=QImage(filename);
}

