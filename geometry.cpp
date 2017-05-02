#include "geometry.h"

int ParticleCount;
bool isReady;

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

    if (index < MAX_NUM_PARTICLES) {
        glPushMatrix();
        glColor4f(0.2, 0.4, 0.6, 0.3);
        glTranslated(position.X() * 2.0, position.Y() * 2.0, position.Z() * 2.0);
        glutSolidSphere(0.02, 10, 10);
        glPopMatrix();
    }

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

ShipBlock::ShipBlock() {
    centerPosition = Vector(0, 0, 0.5);
    acceleration = Vector(0, 0, 0);
    velocity = Vector(0, 0, 0);
    length = 0.12;
    width = 0.09;
    height = 0.1;

    faces.push_back(Plane(Vector(1, 0, 0), Vector(length/2.0, 0, 0)));
    faces.push_back(Plane(Vector(-1, 0, 0), Vector(-length/2.0, 0, 0)));
    faces.push_back(Plane(Vector(0, 1, 0), Vector(0, width/2.0, 0)));
    faces.push_back(Plane(Vector(0, -1, 0), Vector(0, -width/2.0, 0)));
    faces.push_back(Plane(Vector(0, 0, 1), Vector(0, 0, height/2.0)));
    faces.push_back(Plane(Vector(0, 0, -1), Vector(0, 0, -height/2.0)));
}












