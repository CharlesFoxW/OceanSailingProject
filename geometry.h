#ifndef OCEANSAILING_GEOMETRY_H
#define OCEANSAILING_GEOMETRY_H

#include "openGLHeaders.h"
#include <vector>
#include <math.h>

using namespace std;

#define PI 3.141592653589793238462643383279


class Vector {
public:
    Vector() {x = 0; y = 0; z = 0;}
    Vector(double a, double b, double c) {x = a; y = b; z = c;}
    Vector(const Vector& vec) {*this = vec;}    // copy constructor
    // operators to simplify vector calculations:
    Vector& operator=(double scalr) {x = scalr; y = scalr; z = scalr; return *this;}
    Vector& operator+=(const Vector& vec)     { x += vec.x; y += vec.y; z += vec.z; return *this; }
    Vector& operator-=(const Vector& vec)     { x -= vec.x; y -= vec.y; z -= vec.z; return *this; }
    Vector& operator/=(const Vector& vec)     { x /= vec.x; y /= vec.y; z /= vec.z; return *this; }
    Vector& operator/=(const double b)     { x /= b; y /= b; z /= b; return *this; }
    Vector& operator*=(const double b)     { x *= b; y *= b; z *= b;   return *this; }
    Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
    Vector operator-(const Vector &b) const { return Vector(x - b.x,y - b.y,z - b.z); }
    Vector operator*(const double b) const { return Vector(x * b,y * b,z * b); }
    Vector operator/(const double b) const { return Vector(x / b,y / b,z / b); }
    Vector operator*(const Vector& v) const { return Vector(x * v.x, y * v.y, z * v.z); }
    Vector& operator*=(const Vector& v)     { x *= v.x; y *= v.y; z *= v.z; return *this; }

    double X() {return x;}
    double Y() {return y;}
    double Z() {return z;}
    double maxValue() {return x > y && x > z ? x : y > z ? y : z;}
    double dotProduct(Vector vec) {return x * vec.x + y * vec.y + z * vec.z;}
    Vector crossProduct(Vector vec) {return Vector(y * vec.z - vec.y * z, -x * vec.z + vec.x * z, x * vec.y - vec.x * y);}
    double magnitude() {return sqrt(x * x + y * y + z * z);}
    Vector normalize() {
        if (this->magnitude() != 1.0 && this->magnitude() != 0.0) {
            x = x / this->magnitude();
            y = y / this->magnitude();
            z = z / this->magnitude();
        }
        return *this;
    }
    Vector getNormal() {
        Vector vec = *this;
        vec.normalize();
        return vec;
    }

private:
    double x;
    double y;
    double z;
};

class Plane {
public:
    Plane() {}
    Plane(Vector n, Vector p) {normalVector = n.normalize(); point = p;}
    Vector Normal() {return normalVector;}
    Vector Point() {return point;}

private:
    Vector normalVector;
    Vector point;

};


class Particle {
public:
    Particle();
    Particle(Vector p);
    Particle(Vector p, Vector v);

    void draw();
    void clear();

    Vector position;
    Vector velocity;
    Vector force;
    Vector acceleration;
    double density;
    double pressure;
    bool isOnSurface;
    int index;

};

class Grid {
public:
    Grid();
    Grid(int, int, int);
    ~Grid();
    vector<Particle>& operator()(int x, int y, int z) {
        return particlesInGrid[x + y * dimensionX + z * dimensionX * dimensionY];
    }
    vector<Particle>& refParticleVectorAt(int x, int y, int z);

    int dimensionX;
    int dimensionY;
    int dimensionZ;
    int cellCount;
    vector<Particle> * particlesInGrid;
};


#endif //OCEANSAILING_GEOMETRY_H
