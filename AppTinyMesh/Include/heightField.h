#ifndef HEIGHTFIELD_H
#define HEIGHTFIELD_H
#include "mesh.h"
#include "qimage.h"

class HeightField{
private:
    QImage img;
    int echelle;

public:
    HeightField(const QString& filename,const int& echelle);
    ~HeightField();

    Mesh generate();


};

#endif // HEIGHTFIELD_H
