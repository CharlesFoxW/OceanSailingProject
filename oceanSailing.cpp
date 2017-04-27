#include "oceanSailing.h"


OceanSailing::OceanSailing() {

    ParticleCount = 0;
    entranceVelocity = Vector(-2.0, 0, -2.0);
    gravity = Vector(0, 0, GRAVITY_ACCELERATION);

    grid = Grid((int) ceil(SIZE_X / CELL_SIZE), (int) ceil(SIZE_Y / CELL_SIZE), (int) ceil(SIZE_Z / CELL_SIZE));
    boundaries.push_back(Plane(Vector(1, 0, 0), Vector(-SIZE_X / 2.0, 0, 0)));
    boundaries.push_back(Plane(Vector(-1, 0, 0), Vector(SIZE_X / 2.0, 0, 0)));
    boundaries.push_back(Plane(Vector(0, 1, 0), Vector(0, -SIZE_Y / 2.0, 0)));
    boundaries.push_back(Plane(Vector(0, -1, 0), Vector(0, SIZE_Y / 2.0, 0)));
    boundaries.push_back(Plane(Vector(0, 0, 1), Vector(0, 0, -0.1)));

    alterGridByParticle();
}

void OceanSailing::addParticle(Vector p, Vector v) {
    vector<Particle>& firstCell = grid.refParticleVectorAt(0, 0, 0);
    firstCell.push_back(Particle(p, v));
}

void OceanSailing::alterGridByParticle() {

    for (int i = 0; i < grid.dimensionX; i++) {
        for (int j = 0; j < grid.dimensionY; j++) {
            for (int k = 0; k < grid.dimensionZ; k++) {

                vector<Particle>& particleVector = grid.refParticleVectorAt(i, j, k);

                for (int l = 0; l < particleVector.size(); l++) {
                    //Particle ptcl = particleVector[l];

                    int new_x = (int) floor((particleVector[l].position.X() + SIZE_X / 2.0) / CELL_SIZE);
                    int new_y = (int) floor((particleVector[l].position.Y() + SIZE_Y / 2.0) / CELL_SIZE);
                    int new_z = (int) floor((particleVector[l].position.Z() + SIZE_Z / 2.0) / CELL_SIZE);

                    // bound in the cell in grid:
                    if (new_x < 0)
                        new_x = 0;
                    else if (new_x >= grid.dimensionX)
                        new_x = grid.dimensionX - 1;

                    if (new_y < 0)
                        new_y = 0;
                    else if (new_y >= grid.dimensionY)
                        new_y = grid.dimensionY - 1;

                    if (new_z < 0)
                        new_z = 0;
                    else if (new_z >= grid.dimensionZ)
                        new_z = grid.dimensionZ - 1;

                    // If the particle has moved:
                    if (i != new_x || j != new_y || k != new_z) {
                        grid.refParticleVectorAt(new_x, new_y, new_z).push_back(particleVector[l]);
                        particleVector[l] = particleVector.back();
                        particleVector.pop_back();
                        l--;    // re-consider this particle
                    }
                }

            }
        }
    }

}

// computer acceleration for the particle system:
double OceanSailing::kernelPoly6(double radius) {
    return 315.0 / (64.0 * PI * pow(CELL_SIZE, 9)) * pow(CELL_SIZE * CELL_SIZE - radius * radius, 3);
}

Vector OceanSailing::kernelPoly6Gradient(double radius, Vector v) {
    double coeff =  -945.0 / (32.0 * PI * pow(CELL_SIZE, 9))
                    * pow(CELL_SIZE * CELL_SIZE - radius * radius, 2);
    //printf("%.4f\n", coeff);
    return v * coeff;
}

Vector OceanSailing::kernelSpikyGradient(double radius, Vector v) {
    double coeff = -45.0 / (PI * pow(CELL_SIZE, 6)) * pow(CELL_SIZE - radius, 2) / radius;
    return v * coeff;
}

double OceanSailing::kernelPoly6Laplacian(double radius) {
    double coeff = -945.0 / (32.0 * PI * pow(CELL_SIZE, 9));
    return coeff * (CELL_SIZE * CELL_SIZE - radius * radius)
           * (3.0 * CELL_SIZE * CELL_SIZE - 7.0 * radius * radius);
}

double OceanSailing::kernelVicosityLaplacian(double radius) {
    return 45.0 / (PI * pow(CELL_SIZE, 6)) * (CELL_SIZE - radius);
}

void OceanSailing::computerDensityAndPressure() {

    for (int i = 0; i < grid.dimensionX; i++) {
        for (int j = 0; j < grid.dimensionY; j++) {
            for (int k = 0; k < grid.dimensionZ; k++) {

                vector<Particle> &particleVector = grid.refParticleVectorAt(i, j, k);

                for (int l = 0; l < particleVector.size(); l++) {
                    //Particle& ptcl = particleVector[l];
                    particleVector[l].density = 0;  // re-calculate the density each time.
                    // density at particle positions are determined by
                    // the neighbour grids' particles:
                    for (int x = -1; x <= 1; x++) {
                        // don't be outside of the grid:
                        if (i + x < 0 || i + x >= grid.dimensionX)
                            continue;

                        for (int y = -1; y <= 1; y++) {
                            if (j + y < 0 || j + y >= grid.dimensionY)
                                continue;

                            for (int z = -1; z <= 1; z++) {
                                if (k + z < 0 || k + z >= grid.dimensionZ)
                                    continue;

                                vector<Particle>& neighbourParticleVector = grid.refParticleVectorAt(i + x, j + y, k + z);

                                for (int m = 0; m < neighbourParticleVector.size(); m++) {
                                    Particle neighbourPtcl = neighbourParticleVector[m];
                                    Vector distanceVector = particleVector[l].position - neighbourPtcl.position;
                                    double distance = distanceVector.magnitude();
                                    if (distance <= CELL_SIZE){
                                        //printf("%.4f\n", kernelPoly6(distance));
                                        particleVector[l].density += MASS * kernelPoly6(distance);
                                    }
                                }
                            }
                        }
                    }

                    particleVector[l].pressure = AIR_STIFFNESS * (particleVector[l].density - WATER_DENSITY);
                }
            }
        }
    }
}

void OceanSailing::computeAcceleration() {

    computerDensityAndPressure();

    for (int i = 0; i < grid.dimensionX; i++) {
        for (int j = 0; j < grid.dimensionY; j++) {
            for (int k = 0; k < grid.dimensionZ; k++) {

                vector<Particle> &particleVector = grid.refParticleVectorAt(i, j, k);

                for (int l = 0; l < particleVector.size(); l++) {
                    //Particle ptcl = particleVector[l];

                    Vector f_viscosity, f_pressure, f_surface, f_gravity;
                    // parameters of the space field containing water:
                    Vector waterFieldNormal;
                    double waterFieldLaplacian = 0;

                    f_gravity = gravity * particleVector[l].density;

                    // water field on each particle determining by neighbours:
                    for (int x = -1; x <= 1; x++) {
                        // don't be outside of the grid:
                        if (i + x < 0 || i + x >= grid.dimensionX)
                            continue;

                        for (int y = -1; y <= 1; y++) {
                            if (j + y < 0 || j + y >= grid.dimensionY)
                                continue;

                            for (int z = -1; z <= 1; z++) {
                                if (k + z < 0 || k + z >= grid.dimensionZ)
                                    continue;

                                vector<Particle>& neighbourParticleVector = grid.refParticleVectorAt(i + x, j + y, k + z);
                                //printf("%d,\n", (int) neighbourParticleVector.size());
                                for (int m = 0; m < neighbourParticleVector.size(); m++) {
                                    if (particleVector[l].index == neighbourParticleVector[m].index)
                                        continue;

                                    Vector distanceVector = particleVector[l].position - neighbourParticleVector[m].position;
                                    double distance = distanceVector.magnitude();

                                    if (distance <= CELL_SIZE){
                                        Vector poly6Gradient = kernelPoly6Gradient(distance, distanceVector);
                                        Vector spikyGradient = kernelSpikyGradient(distance, distanceVector);
                                        //printf("%.4f, %.4f, %.4f.\n",poly6Gradient.X(), poly6Gradient.Y(), poly6Gradient.Z());
                                        // This algorithm involves particle themselves as neighbours.
                                        // Avoid them when calculating pressure and viscosity.
                                        f_pressure += spikyGradient * (particleVector[l].pressure
                                                                       / pow(particleVector[l].density,2)
                                                                       + neighbourParticleVector[m].pressure
                                                                         / pow(neighbourParticleVector[m].density, 2));
                                        f_viscosity += (neighbourParticleVector[m].velocity - particleVector[l].velocity)
                                                       * kernelVicosityLaplacian(distance) / neighbourParticleVector[m].density;

                                        waterFieldNormal += poly6Gradient / neighbourParticleVector[m].density * MASS;
                                        waterFieldLaplacian += kernelPoly6Laplacian(distance)
                                                               / neighbourParticleVector[m].density * MASS;
                                    }
                                }
                            }
                        }
                    }

                    f_pressure *= -MASS * particleVector[l].density;
                    f_viscosity *= VISCOSITY * MASS;
                    //printf("%.4f, %.4f, %.4f.\n", f_pressure.X(), f_pressure.Y(), f_pressure.Z());
                    //printf("%.4f\n", waterFieldNormal.magnitude());

                    if (waterFieldNormal.magnitude() >= SURFACE_THRESHOLD) {    // the particle is a surface particle
                        particleVector[l].force = waterFieldNormal;
                        particleVector[l].isOnSurface = true;
                        f_surface =  waterFieldNormal.normalize() * waterFieldLaplacian * (-WATER_SURFACE_TENSION);
                        //printf("%.4f, %.4f, %.4f.\n",f_surface.X(), f_surface.Y(), f_surface.Z());
                    } else {
                        particleVector[l].isOnSurface = false;
                    }

                    // from fluid dynamics:
                    particleVector[l].acceleration = (f_pressure + f_viscosity + f_surface + f_gravity) / particleVector[l].density;
                    particleVector[l].acceleration += accFromCollisionDetection(particleVector[l]);
                }
            }
        }
    }

}

Vector OceanSailing::accFromCollisionDetection(Particle p) {
    Vector acc = Vector(0, 0, 0);
    // Using spring system model here:
    for (int i = 0; i < boundaries.size(); i++) {
        Plane& boundary = boundaries[i];
        double distanceOverBoundary = (boundary.Point() - p.position).dotProduct(boundary.Normal()) + 0.005;

        if (distanceOverBoundary > 0) { // outside boundary, apply spring back force
            acc += boundary.Normal() * distanceOverBoundary * BOUNDARY_STIFFNESS;
            acc += boundary.Normal() * p.velocity.dotProduct(boundary.Normal()) * DAMPING_COEFF;
        }
    }
    return acc;
}

void OceanSailing::alterSceneByForwardEuler(double dt) {

    computeAcceleration();

    for (int i = 0; i < grid.cellCount; i++) {

        vector<Particle>& particleVector = grid.particlesInGrid[i];

        for (int j = 0; j < particleVector.size(); j++) {
            particleVector[j].position += particleVector[j].velocity * dt + particleVector[j].acceleration * dt * dt;
            particleVector[j].velocity += particleVector[j].acceleration * dt;
        }
    }
    //printf("count = %d\n", Particle::Count);

    if (ParticleCount < MAX_NUM_PARTICLES) {
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, 0, SIZE_Z*2.0+CELL_SIZE*0.6), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, 0, SIZE_Z*2.0), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, 0, SIZE_Z*2.0+CELL_SIZE*-0.6), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, CELL_SIZE*0.6, SIZE_Z*2.0+CELL_SIZE*0.3), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, CELL_SIZE*-0.6, SIZE_Z*2.0+CELL_SIZE*0.3), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, CELL_SIZE*0.6, SIZE_Z*2.0+CELL_SIZE*-0.3), entranceVelocity);
        addParticle(Vector(SIZE_X/2.0-CELL_SIZE/2.0, CELL_SIZE*-0.6, SIZE_Z*2.0+CELL_SIZE*-0.3), entranceVelocity);
    }

    alterGridByParticle();
}

void OceanSailing::drawScene() {

    for (int i = 0; i < grid.cellCount; i++) {

        vector<Particle>& particleVector = grid.particlesInGrid[i];

        for (int j = 0; j < particleVector.size(); j++) {
            particleVector[j].draw();
        }
    }

}

























