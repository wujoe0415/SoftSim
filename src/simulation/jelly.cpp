#include "jelly.h"

#include "Eigen/Dense"
#include <iostream>
#include "../util/helper.h"
namespace simulation {
constexpr float g_cdK = 2500.0f;
constexpr float g_cdD = 50.0f;

Jelly::Jelly()
    : particleNumPerEdge(10), 
      jellyLength(2.0),
      initialPosition(Eigen::Vector3f(0.0, 0.0, 0.0)),
      springCoefStruct(g_cdK),
      springCoefShear(g_cdK),
      springCoefBending(g_cdK),
      damperCoefStruct(g_cdD),
      damperCoefShear(g_cdD),
      damperCoefBending(g_cdD) {
    particleNumPerFace = particleNumPerEdge * particleNumPerEdge;
    initializeParticle();
    initializeSpring();
}

Jelly::Jelly(const Eigen::Vector3f &a_kInitPos, const float jellyLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      jellyLength(jellyLength),
      initialPosition(a_kInitPos),
      springCoefStruct(dSpringCoef),
      springCoefShear(dSpringCoef),
      springCoefBending(dSpringCoef),
      damperCoefStruct(dDamperCoef),
      damperCoefShear(dDamperCoef),
      damperCoefBending(dDamperCoef) {
    particleNumPerFace = numAtEdge * numAtEdge;
    initializeParticle();
    initializeSpring();
}

int Jelly::getParticleNum() const { return static_cast<int>(particles.size()); }

int Jelly::getSpringNum() const { return static_cast<int>(springs.size()); }

int Jelly::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Jelly::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
    int r = -1;

    switch (a_ciSide) {
        case 1:  // [a_ciI][a_ciJ][0] bottom face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ;
            break;
        case 6:  // [a_ciI][a_ciJ][9] top face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ + particleNumPerEdge - 1;
            break;
        case 2:  // [a_ciI][0][a_ciJ] front face
            r = particleNumPerFace * a_ciI + a_ciJ;
            break;
        case 5:  // [a_ciI][9][a_ciJ] back face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * (particleNumPerEdge - 1) + a_ciJ;
            break;
        case 3:  // [0][a_ciI][a_ciJ] left face
            r = particleNumPerEdge * a_ciI + a_ciJ;
            break;
        case 4:  // [9][a_ciI][a_ciJ] ra_ciIght face
            r = particleNumPerFace * (particleNumPerEdge - 1) + particleNumPerEdge * a_ciI + a_ciJ;
            break;
    }

    return r;
}

Particle &Jelly::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Jelly::getParticlePointer() { return &particles; }

Spring &Jelly::getSpring(int springIdx) { return springs[springIdx]; }

void Jelly::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        springCoefStruct = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        springCoefShear = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        springCoefBending = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::BENDING);
    }
}

void Jelly::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        damperCoefStruct = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        damperCoefShear = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        damperCoefBending = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::BENDING);
    }
}

void Jelly::resetJelly(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Jelly::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}

void Jelly::computeInternalForce() {
    // TODO#2-3: Compute the internal force (including spring force and damper force) for each spring.
    //   1. Read the start-particle and end-particle index from spring.
    //   2. Use `getPosition()` to get particle i's position and use `getVelocity()` get particle i's velocity.
    //   3. Call `computeSpringForce` and `computeDamperForce` to compute spring force and damper force.
    //   4. Compute net internal force and call `addForce` to apply the force onto particles.
    // Note:
    //   1. Direction of the force.

    for (int i = 0; i < getSpringNum();i++) {
        Spring spring = getSpring(i);
        Eigen::Vector3f startPosition = particles[spring.getSpringStartID()].getPosition();
        Eigen::Vector3f endPosition = particles[spring.getSpringEndID()].getPosition();
        Eigen::Vector3f startVelocity = particles[spring.getSpringStartID()].getVelocity();
        Eigen::Vector3f endVelocity = particles[spring.getSpringEndID()].getVelocity();

        Eigen::Vector3f springForce =
            computeSpringForce(startPosition, endPosition, spring.getSpringCoef(), spring.getSpringRestLength());
        Eigen::Vector3f damperForce =
            computeDamperForce(startPosition, endPosition, startVelocity, endVelocity, spring.getDamperCoef());

        particles[spring.getSpringStartID()].addForce(springForce + damperForce);
        particles[spring.getSpringEndID()].addForce(-springForce - damperForce);
    }
}

Eigen::Vector3f Jelly::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    // TODO#2-1: Compute spring force given the two positions of the spring.
    //   1. Review "particles.pptx¡¨ from p.9 - p.13
    //   2. The sample below just set spring force to zero
    Eigen::Vector3f spring_force = Eigen::Vector3f::Zero();
    if ((positionA - positionB).norm() != 0) {
        spring_force = -1 * springCoef * ((positionA - positionB).norm() - restLength) * (positionA - positionB) /
            (positionA - positionB).norm();
    }
    return spring_force;
}

Eigen::Vector3f Jelly::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO#2-2: Compute damper force given the two positions and the two velocities of the spring.
    //   1. Review "particles.pptx¡¨ from p.9 - p.13
    //   2. The sample below just set damper force to zero

    Eigen::Vector3f damper_force = Eigen::Vector3f::Zero();
    if (((positionA - positionB).norm()) != 0) {
        damper_force = -1 * damperCoef * (velocityA - velocityB).dot((positionA - positionB)) *
            (positionA - positionB) /
            ((positionA - positionB).norm() * (positionA - positionB).norm());
    }
    return damper_force;
}

void Jelly::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                particles.push_back(Particle);

            }
        }
    }
}

void Jelly::initializeSpring() {

    // TODO#1: Connect particles with springs.  
    //   1. Consider the type of springs and compute indices of particles which the spring connect to.
    //   2. Compute rest spring length using particle positions.
    //   2. Iterate the particles. Push spring objects into `springs` vector
    // Note:
    //   1. The particles index can be computed in a similar way as below:
    //   ===============================================
    //   0 1 2 3 ... particlesPerEdge
    //   particlesPerEdge + 1 ....
    //   ... ... particlesPerEdge * particlesPerEdge - 1
    //   ===============================================
    // Here is a simple example which connects the structrual springs along z-axis.

    // struct
    // z-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {

                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 1;
                if (k + 1 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                             Spring::SpringType::STRUCT));
                }
                iNeighborID = iParticleID + particleNumPerEdge;
                if (j + 1 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                             Spring::SpringType::STRUCT));
                }
                iNeighborID = iParticleID + particleNumPerFace;
                if (i + 1 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                             Spring::SpringType::STRUCT));
                }
            }
        }
    }
    int bend = 0;
    // bending
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                // z-direction
                int iNeighborID = iParticleID + 2;
                if (k+2 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                             Spring::SpringType::BENDING));
                }
                // y-direction
                iNeighborID = iParticleID + 2 * particleNumPerEdge;
                if (j+2 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                             Spring::SpringType::BENDING));
                }
                // x-direction
                iNeighborID = iParticleID + 2 * particleNumPerFace;
                if (i+2 < particleNumPerEdge) {
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                             Spring::SpringType::BENDING));
                }
            }
        }
    }

    // shear
    int shear = 0;
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                
                for (int inine = -1; inine < 0; inine++) {
                    if (i + inine < 0 || i + inine >= particleNumPerEdge) continue;
                    for (int jnine = -1; jnine < 2; jnine++) {
                        if (j + jnine < 0 || j + jnine >= particleNumPerEdge) continue;
                        for (int knine = -1; knine < 2; knine++) {
                            if (k + knine < 0 || k + knine >= particleNumPerEdge) continue;
                            if (knine == 0 && jnine == 0) continue;

                            int iNeighborID =
                                (i + inine) * particleNumPerFace + (j + jnine) * particleNumPerEdge + (k + knine);
                            Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                            Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                            Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                            float absLength =
                                sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                            springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                     damperCoefShear, Spring::SpringType::SHEAR));
                        }
                    }
                }
                if (j - 1 >= 0 && k-1>=0) {
                    int iNeighborID =
                        i * particleNumPerFace + (j - 1) * particleNumPerEdge + (k - 1);
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = Length.norm();

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                        damperCoefShear, Spring::SpringType::SHEAR));
                }
                if (j - 1 >= 0 && k + 1 < particleNumPerEdge) {
                    int iNeighborID =
                        i * particleNumPerFace + (j - 1) * particleNumPerEdge + (k + 1);
                    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                    float absLength = Length.norm();

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                        damperCoefShear, Spring::SpringType::SHEAR));
                }
            }
        }
    }
}

void Jelly::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Jelly::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
