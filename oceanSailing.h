#ifndef OCEANSAILING_OCEANSAILING_H
#define OCEANSAILING_OCEANSAILING_H

#include "geometry.h"

#define GRAVITY_ACCELERATION    -9.80665
#define CELL_SIZE    0.05
#define SIZE_X  0.5
#define SIZE_Y  0.5
#define SIZE_Z  0.25
#define MASS    0.02
#define AIR_STIFFNESS   3.0
#define WATER_DENSITY   998.3
#define VISCOSITY   3.5
#define SURFACE_THRESHOLD  9.0
#define WATER_SURFACE_TENSION 0.0728
#define BOUNDARY_STIFFNESS  10000.0
#define DAMPING_COEFF   -10.0

extern int ParticleCount;
extern bool isReady;
extern int waveClock;

class OceanSailing {
public:
    OceanSailing();
    ~OceanSailing() {}

    void addParticle(Vector p, Vector v);
    void alterGridByParticle();

    void computerDensityAndPressure();
    void computeAcceleration();
    void alterSceneByForwardEuler(double);
    Vector accFromCollisionDetection(Particle);

    //SPH algorithms from smoothing kernel function:
    double kernelPoly6(double);
    Vector kernelPoly6Gradient(double, Vector);
    Vector kernelSpikyGradient(double, Vector);
    double kernelPoly6Laplacian(double);
    double kernelVicosityLaplacian(double);

    void drawScene();

    Grid grid;
    vector<Plane> boundaries;

    Vector entranceVelocity;
    Vector gravity;
    bool hasShipBlock;
    int numOfTopParticles;
    double currentWaterLevel;
    ShipBlock myShipBlock;
    bool hasWave;
};









#endif //OCEANSAILING_OCEANSAILING_H
