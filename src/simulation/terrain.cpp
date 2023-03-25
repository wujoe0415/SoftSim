#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>
#include <cmath>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    Eigen::Vector3f hole = PlaneTerrain::hole_position;
    float radius = PlaneTerrain::hole_radius;
    for (int i = 0; i < jelly.getParticleNum(); i++)
    {
        auto part = &(jelly.getParticle(i));
        
        if (Eigen::Vector3f::UnitY().dot(part->getPosition()) < eEPSILON && Eigen::Vector3f::UnitY().dot(part->getVelocity()) < 0) {
            if ((part->getPosition() - hole).norm() <= radius)
                continue;
            auto initv = part->getVelocity();
            initv.y() = initv.y() * (-coefResist);
            part->setVelocity(initv);
            }
        if (abs(Eigen::Vector3f::UnitY().dot(part->getPosition())) < eEPSILON && abs(Eigen::Vector3f::UnitY().dot(part->getVelocity())) < eEPSILON) {
            if ((part->getPosition() - hole).norm() <= radius)
                continue;
            Eigen::Vector3f force = part->getForce();
            Eigen::Vector3f contact = -(Eigen::Vector3f::UnitY().dot(force) * Eigen::Vector3f::UnitY());
            if (Eigen::Vector3f::UnitY().dot(force) > 0)
                contact = Eigen::Vector3f::Zero();
            auto vt = part->getVelocity();
            vt.y() = 0;

            Eigen::Vector3f friction = (-coefFriction * (-1 * Eigen::Vector3f::UnitY().dot(force) * vt));
            if (Eigen::Vector3f::UnitY().dot(force) >= 0)
                friction = Eigen::Vector3f::Zero();
            part->addForce(contact + friction);
        }
    }
    // jelly.computeInternalForce();
    
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
}
// BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.8f;
    float radius = BowlTerrain::radius;
    Eigen::Vector3f position = BowlTerrain::position;
    float m2 = BowlTerrain::mass;
    for (int i = 0; i < jelly.getParticleNum(); i++)
    {
        auto part = &(jelly.getParticle(i));
        float m1 = part->getMass();
        if (part->getPosition().y() < 0 && (part->getPosition().x() - position.x()) * (part->getPosition().x() - position.x()) + (part->getPosition().z() - position.z()) * (part->getPosition().z() - position.z()) <= radius * radius)
        {
            Eigen::Vector3f N = (position - part->getPosition()).normalized();
            Eigen::Vector3f p = position - radius * N;
            Eigen::Vector3f vt, vn;
            Eigen::Vector3f originV = part->getVelocity();
            vn = part->getVelocity().dot(-N) * (-N);
            vt = part->getVelocity() - vn;
            if (N.dot((part->getPosition() - p)) < eEPSILON && N.dot(originV) < 0) {
                vn = -coefResist * vn;
                part->setVelocity(vn + vt);
            }
            if (abs(N.dot(part->getPosition() - p)) < eEPSILON && abs(N.dot(originV)) < eEPSILON /* && vt.norm() < 0.1*/) {
                Eigen::Vector3f force = part->getForce();
                Eigen::Vector3f contact = -(N.dot(force) * N);
                if (N.dot(force) > 0)
                    contact = Eigen::Vector3f::Zero();
                Eigen::Vector3f friction = (-coefFriction *
                    (-1 * N.dot(force) * vt));
                if (N.dot(force) >= 0)
                    friction = Eigen::Vector3f::Zero();
                part->addForce(contact + friction);
            }
        }
    }
    
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data. 
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

   
}
}  // namespace simulation
