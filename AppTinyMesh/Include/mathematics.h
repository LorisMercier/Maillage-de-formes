#pragma once

#include <math.h>
#include <ostream>
#include <iostream>

class Math
{
public:
  static constexpr double Clamp(double, double = 0.0, double = 1.0);

  // Minimum and maximum
  static constexpr double Min(double, double);
  static constexpr double Max(double, double);
  static constexpr double Min(double, double, double);
  static constexpr double Max(double, double, double);

  static constexpr double DegreeToRadian(double);
  static constexpr double RadianToDegree(double);
  static constexpr double Pow2(double);
};

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline constexpr double Math::Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Min(double a, double b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Max(double a, double b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Max(double a, double b, double c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Min(double a, double b, double c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline constexpr double Math::DegreeToRadian(double a)
{
  return a * 3.14159265358979323846 / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline constexpr double Math::RadianToDegree(double a)
{
  return a * 180.0 / 3.14159265358979323846;
}

inline constexpr double Math::Pow2(double a){
    return a*a;
}

// Class
class Vector
{
protected:
  double c[3]; //!< Components.
public:
  //! Empty 
  Vector() {}

  explicit Vector(double);
  explicit Vector(double, double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vector operator+ () const;
  Vector operator- () const;

  // Assignment operators
  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (const Vector&);
  Vector& operator/= (const Vector&);
  Vector& operator*= (double);
  Vector& operator/= (double);

  // Binary operators
  friend int operator> (const Vector&, const Vector&);
  friend int operator< (const Vector&, const Vector&);

  friend int operator>= (const Vector&, const Vector&);
  friend int operator<= (const Vector&, const Vector&);

  // Binary operators
  friend Vector operator+ (const Vector&, const Vector&);
  friend Vector operator- (const Vector&, const Vector&);

  friend constexpr double operator* (const Vector&, const Vector&);

  friend Vector operator* (const Vector&, double);
  friend Vector operator* (double, const Vector&);
  friend Vector operator/ (const Vector&, double);

  friend Vector operator/ (const Vector&, const Vector&);

  // Boolean functions
  friend int operator==(const Vector&, const Vector&);
  friend int operator!=(const Vector&, const Vector&);

  // Norm
  friend double Norm(const Vector&);
  friend double SquaredNorm(const Vector&);

  friend void Normalize(Vector&);
  friend Vector Normalized(const Vector&);

  // Compare functions
  static Vector Min(const Vector&, const Vector&);
  static Vector Max(const Vector&, const Vector&);

  // Abs
  friend Vector Abs(const Vector&);

  // Orthogonal and orthonormal vectors
  Vector Orthogonal() const;
  void Orthonormal(Vector&, Vector&) const;

  friend Vector Lerp(const Vector&, const Vector&, double);
  static Vector Bilinear(const Vector&, const Vector&, const Vector&, const Vector&, double, double);

  // Scale
  Vector Scaled(const Vector&) const;
  Vector Inverse() const;

  friend std::ostream& operator<<(std::ostream&, const Vector&);

public:
  static const Vector Null; //!< Null vector.
  static const Vector X; //!< Vector(1,0,0).
  static const Vector Y; //!< Vector(0,1,0).
  static const Vector Z; //!< Vector(0,0,1).
};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vector::Vector(double a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vector::Vector(double a, double b, double c)
{
  Vector::c[0] = a;
  Vector::c[1] = b;
  Vector::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline double& Vector::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline double Vector::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vector Vector::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vector Vector::operator- () const
{
  return Vector(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vector& Vector::operator+= (const Vector& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vector& Vector::operator-= (const Vector& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vector& Vector::operator*= (double a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Scale a vector.
\param a Scaling vector.
*/
inline Vector Vector::Scaled(const Vector& a) const
{
  return Vector(c[0] * a[0], c[1] * a[1], c[2] * a[2]);
}

/*!
\brief Inverse of a vector.

This function inverses the components of the vector. This is the same as:
\code
Vector v=Vector(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vector Vector::Inverse() const
{
  return Vector(1.0 / c[0], 1.0 / c[1], 1.0 / c[2]);
}

//! Destructive division by a scalar.
inline Vector& Vector::operator/= (double a)
{
  c[0] /= a; c[1] /= a; c[2] /= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.

This is the same as Scale:
\code
Vector u(2.0,-1.0,1.0);
u=u.Scaled(Vector(3.0,1.0,2.0)); // u*=Vector(3.0,1.0,2.0);
\endcode
*/
inline Vector& Vector::operator*= (const Vector& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Destructively divide the components of a vector by another vector.
inline Vector& Vector::operator/= (const Vector& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1]; c[2] /= u.c[2];
  return *this;
}

//! Compare two vectors.
inline int operator> (const Vector& u, const Vector& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]) && (u.c[2] > v.c[2]));
}

//! Compare two vectors.
inline int operator< (const Vector& u, const Vector& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]) && (u.c[2] < v.c[2]));
}

//! Overloaded
inline int operator>= (const Vector& u, const Vector& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]) && (u.c[2] >= v.c[2]));
}

//! Overloaded
inline int operator<= (const Vector& u, const Vector& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]) && (u.c[2] <= v.c[2]));
}

//! Adds up two vectors.
inline Vector operator+ (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vector operator- (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

//! Scalar product.
inline constexpr double operator* (const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}

//! Right multiply by a scalar.
inline Vector operator* (const Vector& u, double a)
{
  return Vector(u.c[0] * a, u.c[1] * a, u.c[2] * a);
}

//! Left multiply by a scalar.
inline Vector operator* (double a, const Vector& v)
{
  return v * a;
}

//! Cross product.
inline Vector operator/ (const Vector& u, const Vector& v)
{
  return Vector(u.c[1] * v.c[2] - u.c[2] * v.c[1], u.c[2] * v.c[0] - u.c[0] * v.c[2], u.c[0] * v.c[1] - u.c[1] * v.c[0]);
}

//! Left multiply by a scalar
inline Vector operator/ (const Vector& u, double a)
{
  return Vector(u.c[0] / a, u.c[1] / a, u.c[2] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vector& u, const Vector& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline int operator!= (const Vector& u, const Vector& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a vector.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a vector instead.
\param u %Vector.
\sa SquaredNorm
*/
inline double Norm(const Vector& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vector.
\sa Norm
*/
inline double SquaredNorm(const Vector& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Return a normalized vector.

Compute the inverse of its norm and scale the components.

This function does not check if the vector is null.
\param u %Vector.
*/
inline Vector Normalized(const Vector& u)
{
  return u * (1.0 / Norm(u));
}

/*!
\brief Computes the absolute value of a vector.
\param u %Vector.
*/
inline Vector Abs(const Vector& u)
{
  return Vector(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

/*!
\brief Return a vector with coordinates set to the minimum coordinates
of the two argument vectors.
*/
inline Vector Vector::Min(const Vector& a, const Vector& b)
{
  return Vector(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1], a[2] < b[2] ? a[2] : b[2]);
}

/*!
\brief Return a vector with coordinates set to the maximum coordinates
of the two argument vectors.
*/
inline Vector Vector::Max(const Vector& a, const Vector& b)
{
  return Vector(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1], a[2] > b[2] ? a[2] : b[2]);
}

/*!
\brief Linear interpolation between two vectors.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vector Lerp(const Vector& a, const Vector& b, double t)
{
  return a + t * (b - a);
}

/*!
\brief Bi-linear interpolation between four vectors.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated vectors.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vector Vector::Bilinear(const Vector& a00, const Vector& a10, const Vector& a11, const Vector& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}




class Matrix{
private:
    double m[9];
public:

    Matrix(); //OK
    Matrix(const Matrix&); //OK
    Matrix(double,double,double,double,double,double,double,double,double); //OK

    ~Matrix(); //OK

    // Access members
    double& operator[] (int); //OK
    double operator[] (int) const; //OK


    Matrix& operator= (const Matrix&);
    Matrix& operator+= (const Matrix&);
    Matrix& operator-= (const Matrix&);
    Matrix& operator*= (const Matrix&);
    Matrix& operator*= (double);
    Matrix& operator/= (double);


    Matrix operator+ (const Matrix&) const;
    Matrix operator- (const Matrix&) const;
    Matrix operator* (const Matrix&) const;
    Matrix operator* (double) const;
    Vector operator* (const Vector&) const;
    Matrix operator/ (double) const;

    Matrix identite();
    //void transpose();
    Matrix transpose();
    Matrix inverse();

    static const Matrix MatRotX(const double&);
    static const Matrix MatRotZ(const double&);
    static const Matrix MatRotY(const double&);

    static const Matrix MatHomothetie(const double&,const double&,const double&);

    friend std::ostream& operator<<(std::ostream&, const Vector&);
};

inline std::ostream& operator<<(std::ostream& os, const Matrix& _m){
    int a = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            os << _m[a] << " ";
            a++;
        }
        os << std::endl;
    }
    return os;
}

/*!
\brief Constructeur d'une matrice égal à 0
*/
inline Matrix::Matrix(){
    for(int i=0;i<9;i++){
        m[i]=0;
    }
}

/*!
\brief Constructeur par copie d'une matrice
\param _m Matrix.
*/
inline Matrix::Matrix(const Matrix & _m){
    for(int i=0;i<8;i++){
        m[i]=_m[i];
    }
}

/*!
 * \brief Constructeur avec valeurs passées en paramètres
 * \param i1 double
 * \param i2 double
 * \param i3 double
 * \param i4 double
 * \param i5 double
 * \param i6 double
 * \param i7 double
 * \param i8 double
 * \param i9 double
 */
inline Matrix::Matrix(double i1,double i2,double i3,double i4,double i5,double i6,double i7,double i8,double i9){
  m[0]=i1;
  m[1]=i2;
  m[2]=i3;
  m[3]=i4;
  m[4]=i5;
  m[5]=i6;
  m[6]=i7;
  m[7]=i8;
  m[8]=i9;
}

inline Matrix::~Matrix(){}

/*!
 * \brief Surcharge de l'opérateur []
 * \param i indice de la valeur à obtenir
 * \return valeur : int
 */
inline double Matrix::operator[] (int i) const {
  return m[i];
}

/*!
 * \brief Surcharge de l'opérateur []
 * \param i indice de la valeur à obtenir
 * \return valeur : int
 */
inline double& Matrix::operator[] (int i){
  return m[i];
}

/*!
 * \brief Surcharge de l'opérateur d'affectation = pour affecter une matrice
 * \param _m Matrix
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator=(const Matrix &_m){
    if(this == &_m ){
        return *this;
    }

    for(int i=0;i<9;i++){
        m[i]=_m[i];
    }

    return *this;
}

/*!
 * \brief Surcharge de l'opérateur += pour additionner deux matrices
 * \param _m Matrix
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator+=(const Matrix &_m){
  for(int i=0; i<9; i++){
    m[i] = m[i] + _m[i];
  }
  return *this;
}

/*!
 * \brief Surcharge de l'opérateur -= pour soustraire deux matrices
 * \param _m Matrix
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator-=(const Matrix &_m){
  for(int i=0; i<9; i++){
    m[i] = m[i] - _m[i];
  }
  return *this;
}

/*!
 * \brief Surcharge de l'opérateur *= pour multipler deux matrices
 * \param _m Matrix
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator*=(const Matrix &_m){
  Matrix mreturn;
  mreturn[0]  = m[0]*_m[0] + m[1]*m[3] + m[2]*m[6];
  mreturn[1]  = m[0]*_m[1] + m[1]*m[4] + m[2]*m[7];
  mreturn[2]  = m[0]*_m[2] + m[1]*m[5] + m[2]*m[8];

  mreturn[3]  = m[3]*_m[0] + m[4]*m[3] + m[5]*m[6];
  mreturn[4]  = m[3]*_m[1] + m[4]*m[4] + m[5]*m[7];
  mreturn[5]  = m[3]*_m[2] + m[4]*m[5] + m[5]*m[8];

  mreturn[6]  = m[6]*_m[0] + m[7]*m[3] + m[8]*m[6];
  mreturn[7]  = m[6]*_m[1] + m[7]*m[4] + m[8]*m[7];
  mreturn[8]  = m[6]*_m[2] + m[7]*m[5] + m[8]*m[8];

  *this = mreturn;

  return *this;
}

/*!
 * \brief Surcharge de l'opérateur *= pour multiplier une matrice à un double
 * \param d double
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator*=(double d){
    for(int i=0; i<9; i++){
        m[i] = m[i] * d;
    }
    return *this;
}

/*!
 * \brief Surcharge de l'opérateur /= pour diviser une matrice par un double
 * \param d double
 * \return Matrix actuelle modifiée
 */
inline Matrix& Matrix::operator/=(double d){
    for(int i=0; i<9; i++){
        m[i] = m[i] / d;
    }
    return *this;
}

/*!
 * \brief Surcharge de l'opérateur + pour additionner deux matrices
 * \param _m Matrix
 * \return nouvelle Matrix
 */
inline Matrix Matrix::operator+(const Matrix &_m) const{
  Matrix mreturn;
  for(int i=0; i<9; i++){
    mreturn[i] = m[i] + _m[i];
  }
  return mreturn;
}

inline Vector Matrix::operator* (const Vector& u) const
{
  return Vector(u[0] *m[0]+u[1] *m[1]+u[2] *m[2],u[0] *m[3]+u[1] *m[4]+u[2] *m[5], u[0] *m[6]+u[1] *m[7]+u[2] *m[8]);
}

/*!
 * \brief Surcharge de l'opérateur - pour soustraire deux matrices
 * \param _m Matrix
 * \return nouvelle Matrix
 */
inline Matrix Matrix::operator-(const Matrix &_m) const{
  Matrix mreturn;
  for(int i=0; i<9; i++){
    mreturn[i] = m[i] - _m[i];
  }
  return mreturn;
}

/*!
 * \brief Surcharge de l'opérateur * pour multiplier deux matrices
 * \param _m Matrix
 * \return nouvelle Matrix
 */
inline Matrix Matrix::operator*(const Matrix &_m) const{
  Matrix mreturn;

  mreturn[0]  = m[0]*_m[0] + m[1]*m[3] + m[2]*m[6];
  mreturn[1]  = m[0]*_m[1] + m[1]*m[4] + m[2]*m[7];
  mreturn[2]  = m[0]*_m[2] + m[1]*m[5] + m[2]*m[8];

  mreturn[3]  = m[3]*_m[0] + m[4]*m[3] + m[5]*m[6];
  mreturn[4]  = m[3]*_m[1] + m[4]*m[4] + m[5]*m[7];
  mreturn[5]  = m[3]*_m[2] + m[4]*m[5] + m[5]*m[8];

  mreturn[6]  = m[6]*_m[0] + m[7]*m[3] + m[8]*m[6];
  mreturn[7]  = m[6]*_m[1] + m[7]*m[4] + m[8]*m[7];
  mreturn[8]  = m[6]*_m[2] + m[7]*m[5] + m[8]*m[8];

  return mreturn;
}

/*!
 * \brief Surcharge de l'opérateur * pour multiplier une matrice et un double
 * \param _m Matrix
 * \return nouvelle Matrix
 */
inline Matrix Matrix::operator*(double d) const{
  Matrix mreturn;
  for(int i=0; i<9; i++){
    mreturn[i] = m[i] * d;
  }
  return mreturn;
}

/*!
 * \brief Surcharge de l'opérateur / pour diviser une matrice par un double
 * \param _m Matrix
 * \return nouvelle Matrix
 */
inline Matrix Matrix::operator/(double d) const{
  Matrix mreturn;
  for(int i=0; i<9; i++){
    mreturn[i] = m[i] / d;
  }
  return mreturn;
}

/*!
 * \brief Créer une matrice identite
 * \return Matrix
 */
inline Matrix identite(){
    return Matrix(1,0,0,0,1,0,0,0,1);
}

/*
inline void Matrix::transpose(){
  double tmp;
  tmp=m[1];
  m[1]=m[3];
  m[3]=tmp;
  tmp=m[2];
  m[2]=m[6];
  m[6]=tmp;
  tmp=m[5];
  m[5]=m[7];
  m[7]=tmp;
}
*/

/*!
 * \brief Transpose la matrice dans une nouvelle matrice
 * \return Matrix
 */
inline Matrix Matrix::transpose(){
    return Matrix(m[0],m[3],m[6],m[1],m[4],m[7],m[2],m[5],m[8]);
}

/*!
 * \brief Inverse la matrice dans une nouvelle matrice
 * \return Matrix
 */
inline Matrix Matrix::inverse(){
    double c11,c12,c13,c21,c22,c23,c31,c32,c33;
    double det;

    //Calcul des cofacteurs
    c11 = m[4] * m[8] - m[5] * m[7];
    c12 = m[3] * m[8] - m[5] * m[6];
    c13 = m[3] * m[7] - m[4] * m[6];

    c21 = m[1] * m[8] - m[2] * m[7];
    c22 = m[0] * m[8] - m[2] * m[6];
    c23 = m[0] * m[7] - m[1] * m[6];

    c31 = m[1] * m[5] - m[2] * m[4];
    c32 = m[0] * m[5] - m[2] * m[3];
    c33 = m[0] * m[4] - m[1] * m[3];

    //Complémentaire = Transposé de la matrice des cofacteurs
    Matrix comp(c11,-c12,c13,-c21,c22,-c23,c31,-c32,c33);

    comp = comp.transpose();

    //Calcul du déterminant via le développement de cofacteur
    det = m[0]*c11 - m[1] * c12 + m[2] * c13;

    if(det != 0){
        return comp/det;
    }
    else{
        std::cout<<"Erreur inversion : determinant nul";
        return *this;
    }
}

inline const Matrix Matrix::MatRotX(const double& d){
    return Matrix(1,0,0,0,cos(d),-sin(d),0,sin(d),cos(d));
}

 inline const Matrix Matrix::MatRotY(const double& d){
    return Matrix(cos(d),0,-sin(d),0,1,0,sin(d),0,cos(d));
}

 inline const Matrix Matrix::MatRotZ(const double& d){
    return Matrix(cos(d),-sin(d),0,sin(d),cos(d),0,0,0,1);
}

 inline const Matrix Matrix::MatHomothetie(const double& x,const double& y,const double& z){
    return Matrix(x,0,0,0,y,0,0,0,z);
}

