#include "mesh.h"

/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize 
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}

/*!
 * \brief Creates a disc.
 * \param disque The disc.
 * \param nbpoint Number of points drawing the circumference of the disc.
 */
Mesh::Mesh(const Disque& disque, const int& nbpoint)
{
  // Vertices
  vertices.resize(nbpoint + 1);
  float alpha;
  float step = 2.0 * M_PI / (nbpoint);

  vertices[0] = disque.c;

  // Reserve space for the triangle array
  varray.reserve(nbpoint * 3);
  narray.reserve(nbpoint * 3);

  // Normals
  ///+ disque.c pour régler normal ?
  normals.push_back(Vector(0, 0, 1));

  for (int i = 0; i < nbpoint; i++)
  {
      alpha = i*step;
      vertices[i+1] = Vector(cos(alpha),sin(alpha),0)*disque.r + disque.c;
      AddTriangle(0, 1+(i+1)%nbpoint, i+1, 0);
  }
}

/*!
 * \brief Creates a cylinder.
 * \param cylindre The cylinder.
 * \param nbpoint Number of points drawing the circumference of each disk.
 * \param nbCercle Number of circle used.
 * \param afficherDisc Boolean to display or not the disks at the end.
 */
Mesh::Mesh(const Cylindre& cylindre, const int& nbpoint,const int& nbCercle, const bool& afficherDisc)
{
    float alpha;
    float step = 2.0 * M_PI / (nbpoint);
    float steph = cylindre.h/(nbCercle-1);
    Vector c=cylindre.c-Vector(0,0,(cylindre.h/2));


  // Reserve space for the triangle array
  varray.reserve(nbpoint*nbCercle*6);
  narray.reserve(nbpoint*nbCercle*6);

  // Création des disques
  for (int i = 0; i < nbCercle; ++i) {
      for (int j = 0; j < nbpoint; j++)
      {
          alpha = j*step;
          vertices.push_back(Vector(cos(alpha),sin(alpha),0)*cylindre.r + c + Vector(0,0,steph*i));
          normals.push_back(vertices[(j+nbpoint*i)] - (c + Vector(0,0,steph*i)));
      }
  }


  //Gestion des centres des 2 disques
  if(afficherDisc){
      vertices.push_back(c);
      vertices.push_back(c + Vector(0,0,cylindre.h));

      normals.push_back(Vector(0, 0, -1));
      normals.push_back(Vector(0, 0, 1));
  }
  // Triangle de côté
  for (int i = 0; i < nbCercle - 1; ++i) {
      for(int j = 0; j < nbpoint; j++){
          AddSmoothTriangle(j+nbpoint*i,j+nbpoint*i, j+nbpoint*(i+1),j+nbpoint*(i+1), (j+1)%nbpoint+nbpoint*(i+1),(j+1)%nbpoint+nbpoint*(i+1));
          AddSmoothTriangle(j+nbpoint*i,j+nbpoint*i, (j+1)%nbpoint+nbpoint*(i+1),(j+1)%nbpoint+nbpoint*(i+1), (j+1)%nbpoint+nbpoint*i,(j+1)%nbpoint+nbpoint*i);
      }
  }

  // Triangle des disques des extremitées
  if(afficherDisc){
      for(int j = 0; j < nbpoint; j++){
        AddSmoothTriangle(nbpoint*nbCercle, nbpoint*nbCercle,j, nbpoint*nbCercle, (j+1)%nbpoint, nbpoint*nbCercle);
        AddSmoothTriangle(nbpoint*nbCercle+1, nbpoint*nbCercle+1,j+nbpoint*(nbCercle-1), nbpoint*nbCercle+1, (j+1)%nbpoint+nbpoint*(nbCercle-1), nbpoint*nbCercle+1);
    }
  }
}

/*!
 * \brief Creates a sphere.
 * \param sphere The sphere.
 * \param nbPointParCercle Number of points drawing the circumference of each disk.
 * \param nbCercle Number of circle used.
 * \param demiSphere Boolean to create a half-sphere instead of a sphere.
 */
Mesh::Mesh(const Sphere& sphere, const int& nbPointParCercle, const int& nbCercle,const bool& demiSphere)
{
  float alpha,beta;
  int divBeta = nbPointParCercle;


  // Vertices
  vertices.resize((nbCercle) * divBeta + 2);
  normals.resize((nbCercle) * divBeta + 2);

  // Reserve space for the triangle array
  varray.reserve((nbCercle*divBeta) * 6);
  narray.reserve((nbCercle*divBeta) * 6);

  // Création des Vertices sur les cercles intermédiaires
  for (int i = 0; i < nbCercle; i++)
  {
      if(demiSphere)
          alpha = -0.5 * M_PI + double(i+1) * M_PI / (nbCercle)/2 ;
      else
          alpha = -0.5 * M_PI + double(i+1) * M_PI / (nbCercle+1) ;


      for (int j = 0; j < divBeta; j++)
      {
          beta = (double(j) * 2.0 * M_PI) / divBeta;
          vertices[j + i*divBeta] = Vector(cos(alpha)*cos(beta), cos(alpha)*sin(beta), -sin(alpha)) * sphere.r + sphere.c;
          normals[j + i*divBeta] = vertices[j + i*divBeta] - sphere.c;
      }
  }

  // Création des Vertices SOMMET
  vertices[(nbCercle)*divBeta] = Vector(0,0,1) * sphere.r + sphere.c;
  if(demiSphere)
      vertices[(nbCercle)*divBeta + 1] = sphere.c;
  else
      vertices[(nbCercle)*divBeta + 1] = Vector(0,0,-1)* sphere.r + sphere.c;

  normals[(nbCercle)*divBeta] = Vector(0,0,1);
  normals[(nbCercle)*divBeta + 1] = Vector(0,0,-1);

  // Remplissage triangles cercles intermédiaires
  for (int i = 0; i < nbCercle - 1; i++)
  {
      for (int j = 0; j < divBeta; j++)
      {
          AddSmoothTriangle(j + i*divBeta,j + i*divBeta, (j+1)%divBeta + i*divBeta,(j+1)%divBeta + i*divBeta, j+divBeta + i*divBeta, j+divBeta + i*divBeta);
          AddSmoothTriangle((j+1)%divBeta + i*divBeta,(j+1)%divBeta + i*divBeta, (j+1)%divBeta+divBeta + i*divBeta,(j+1)%divBeta+divBeta + i*divBeta, j+divBeta + i*divBeta, j+divBeta + i*divBeta);
      }
  }

  // Remplissage triangles sommet
  for (int j = 0; j < divBeta; j++)
  {
      AddSmoothTriangle((nbCercle)*divBeta,(nbCercle)*divBeta, (j+1)%divBeta, (j+1)%divBeta, j,j);
      if(!demiSphere)
          AddSmoothTriangle((nbCercle)*divBeta + 1, (nbCercle)*divBeta + 1, (j+1)%divBeta + (nbCercle - 1)*divBeta,(j+1)%divBeta + (nbCercle - 1)*divBeta, j + (nbCercle-1) * divBeta,j + (nbCercle-1) * divBeta);
  }
}

/*!
 * \brief Creates a torus.
 * \param tore The torus.
 * \param nbpoint Number of points drawing the circumference of each disk.
 * \param nbCercle Number of circle used.
 */
Mesh::Mesh(const Tore& tore, const int& nbpoint,const int& nbCercle)
{
  float alpha,alpha2;
  float step = 2.0 * M_PI / (nbpoint);
  float step2 = 2.0 * M_PI / (nbCercle);
  Vector c;


  // Reserve space for the triangle array
  varray.reserve(nbpoint*nbCercle*6);
  narray.reserve(nbpoint*nbCercle*6);


    for(int i=0;i<nbCercle;i++){
        alpha2 = i*step2;
        c = Vector(cos(alpha2),sin(alpha2),0)*tore.r + tore.c;
      for (int j = 0; j < nbpoint; j++)
      {
          alpha = j*step;
          vertices.push_back((Matrix::MatRotZ(alpha2)*Vector(cos(alpha),0,sin(alpha)))*tore.rtube+c);
          normals.push_back(vertices.back() - c);
      }
    }

    for (int i = 0; i < nbCercle; i++)
    {
        for (int j = 0; j < nbpoint; ++j) {
            AddSmoothTriangle(i*nbpoint+j,i*nbpoint+j, (j%nbpoint)+ ((i+1)%nbCercle)*nbpoint, (j%nbpoint)+ ((i+1)%nbCercle)*nbpoint, (j+1)%nbpoint+(i*nbpoint),(j+1)%nbpoint+(i*nbpoint));
            AddSmoothTriangle((j%nbpoint)+((i+1)%nbCercle)*nbpoint,(j%nbpoint)+((i+1)%nbCercle)*nbpoint, (j+1)%nbpoint+(((i+1)%nbCercle)*nbpoint), (j+1)%nbpoint+(((i+1)%nbCercle)*nbpoint), (j+1)%nbpoint+(i*nbpoint),(j+1)%nbpoint+(i*nbpoint));
        }

    }
}

/*!
 * \brief Creates a capsule.
 * \param capsule The capsule.
 * \param nbpointCyl Number of points drawing the circumference of each disk.
 * \param nbCercleCyl Number of circles used by the cylinder.
 * \param nbCercleSpher Number of circles used by the half-sphere.
 */
Mesh::Mesh(const Capsule& capsule, const int& nbpoint, const int& nbCercleCyl, const int& nbCercleSpher)
{
    this->Merge(Mesh(Sphere(capsule.c+Vector(0,0,capsule.h/2-capsule.r),capsule.r),nbpoint,nbCercleSpher,true));
    this->Merge(Mesh(Cylindre(capsule.c,capsule.r,capsule.h-capsule.r*2),nbpoint,nbCercleCyl,false));
    Mesh caps=Mesh(Sphere(capsule.c+Vector(0,0,capsule.h/2-capsule.r),capsule.r),nbpoint,nbCercleSpher,true);
    caps.transform(Matrix::MatRotX(M_PI));
    this->Merge(caps);
}


/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}



#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
 * \brief Add vertices and normals to the existing one to merge two meshes.
 * \param b The mesh to add.
 */
void Mesh::Merge(const Mesh& b){

    int sn = normals.size();
    int sv = vertices.size();
    for (int i = 0; i < b.narray.size(); i++) {
        narray.push_back(b.narray[i]+sn);
    }
    for (int i = 0; i < b.varray.size(); i++) {
        varray.push_back(b.varray[i]+sv);
    }
    for (int i = 0; i < b.normals.size(); i++) {
        normals.push_back(b.normals[i]);
    }
    for (int i = 0; i < b.vertices.size(); i++) {
        vertices.push_back(b.vertices[i]);
    }

}

/*!
 * \brief Apply a matrix on all vertices and normals.
 * \param m The matrix to apply.
 */
void Mesh::transform(const Matrix& m){
    for(int i=0;i<vertices.size();i++){
        vertices[i]=(m)*vertices[i];
    }
    for(int i=0;i<normals.size();i++){
        normals[i]=m*normals[i];
    }

}

/*!
 * \brief Apply a translation on each of the points of the mesh.
 * \param v Translation vector.
 */
void Mesh::transfer(const Vector& v){
    for(int i=0;i<vertices.size();i++){
        vertices[i]=v+vertices[i];
    }
}

/*!
 * \brief Apply a deformation to the points in the circle according to the distance from the center of the circle.
 * \param s The sphere.
 * \param t The deformation vector.
 */
void Mesh::SphereWarp(const Sphere& s,const Vector& t){
    double d;
    Vector v;
    for(int i=0;i<vertices.size();i++){
        d=Math::Pow2(s.c[0]-vertices[i][0])+Math::Pow2(s.c[1]-vertices[i][1])+Math::Pow2(s.c[2]-vertices[i][2]);
        if(d<Math::Pow2(s.r)){
            v=t*(1-(d/Math::Pow2(s.r)));
            vertices[i]=vertices[i]+v;
        }
    }

}

/*!
 * \brief Calculates if the ray intersects the mesh.
 * \param r The ray.
 * \return True if the half line intersects otherwise false.
 */
bool Mesh::Intersect(const Ray& r)const{
    double x,y,z;
    for (int i = 0; i < this->Triangles(); ++i) {
        Triangle t=this->GetTriangle(i);
        if(t.Intersect(r,x,y,z)&&x>0){
            return true;
        }
    }
    return false;
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}

