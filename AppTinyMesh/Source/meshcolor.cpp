#include "meshcolor.h"

/*!
\brief Create an empty mesh.
*/
MeshColor::MeshColor()
{
}

/*!
\brief Constructor from a Mesh with color array and indices.
\param m Base mesh.
\param cols Color array.
\param carr Color indexes, should be the same size as Mesh::varray and Mesh::narray.
*/
MeshColor::MeshColor(const Mesh& m, const std::vector<Color>& cols, const std::vector<int>& carr) : Mesh(m), colors(cols), carray(carr)
{
}

/*!
\brief Constructor from a Mesh.
\param m the base mesh
*/
MeshColor::MeshColor(const Mesh& m) : Mesh(m)
{
	colors.resize(vertices.size(), Color(1.0, 1.0, 1.0));
    carray = varray;
}

/*!
\brief Constructor from a Mesh.
\param m the base mesh
*/
MeshColor::MeshColor(const Mesh& m,const Color& c) : Mesh(m)
{
    colors.resize(vertices.size(), c);
    carray = varray;
}

/*!
\brief Empty.
*/
MeshColor::~MeshColor()
{
}

/*!
 * \brief MeshColor::Accessibility calculates the accessibility of a vertex
 * \param indexVertice index of the vertex in the table vertices
 * \param nbrayon number of radius from the vertex that will be intersected or not to the mesh
 * \return double designant le taux d'accessibilité
 */
double MeshColor::Accessibility(const int& indexVertice,const int& nbrayon){
    double x,y,z;
    Vector vv;
    Ray r;
    double res = 0.0;
    Vector v = vertices[indexVertice];
    Vector n = normals[indexVertice];
    for(int i=0;i<nbrayon;i++){
            x = ((rand()%20)/10.0)-1;
            y = ((rand()%20)/10.0)-1;
            z = ((rand()%20)/10.0)-1;
            vv = Vector(x,y,z);
        if(Math::RadianToDegree(acos((vv[0]*n[0]+vv[1]*n[1]+vv[2]*n[2]) / (Norm(n)*Norm(vv))))>90)
            vv*=-1;
        r= Ray(v+0.0001*n,Normalized(vv));
        if(!this->Intersect(r)){
            res+=1.0;
        }
    }
    return res/nbrayon;
}

/*!
 * \brief MeshColor::CalculOcclusion calculates the accessibility of each vertex and changes the associated color
 * \param nbrayon number of radius from the vertex that will be intersected or not to the mesh
 */
void MeshColor::CalculOcclusion(const int& nbrayon){
    double res;
    srand (time(NULL));
    for (int i = 0; i < this->Vertexes(); ++i) {
        res=Accessibility(i,nbrayon);
        colors[i]=colors[i]*res;

    }

}
