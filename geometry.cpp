#include "geometry.h"

int ParticleCount;

Particle::Particle() {
}

Particle::Particle(Vector p) {
    position = p;
    index = ParticleCount++;
}

Particle::Particle(Vector p, Vector v) {
    position = p;
    velocity = v;
    index = ParticleCount++;
}

void Particle::draw() {

    glPushMatrix();
    //if (isOnSurface)
        //glColor4f(0.2, 0.5, 0.4, 0.2);
    //else
        glColor4f(0.2, 0.5, 0.6, 0.3);
    glTranslated(position.X() * 2.0, position.Y() * 2.0, position.Z() * 2.0);
    glutSolidSphere(0.025, 10, 10);
    glPopMatrix();

    if (isOnSurface) {

    }
}

void Particle::clear() {
    position = Vector();
    velocity = Vector();
    force = Vector();
    acceleration = Vector();
    density = 0;
    pressure = 0;
}

Grid::Grid() {
    dimensionX = 0;
    dimensionY = 0;
    dimensionZ = 0;
    cellCount = 0;
}

Grid::Grid(int x, int y, int z) {
    dimensionX = x;
    dimensionY = y;
    dimensionZ = z;
    cellCount = x * y * z;

    particlesInGrid = new vector<Particle>[cellCount];
}

Grid::~Grid() {
}

vector<Particle>& Grid::refParticleVectorAt(int x, int y, int z) {
    return particlesInGrid[x + y * dimensionX + z * dimensionX * dimensionY];
}











