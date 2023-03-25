#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    for (int index = 0; index < particleSystem.jellyCount; index++) {
        auto jelly = particleSystem.getJellyPointer(index);

        for (int i = 0; i < jelly->getParticleNum(); i++) {
            jelly->getParticle(i).addPosition(jelly->getParticle(i).getVelocity() * particleSystem.deltaTime);
            jelly->getParticle(i).addVelocity(jelly->getParticle(i).getAcceleration() * particleSystem.deltaTime);

            jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*jelly);
    }
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.15 - p.16
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    for (int index = 0; index < particleSystem.jellyCount; index++) {
        auto jelly = particleSystem.getJellyPointer(index);
        std::vector<Particle> origin;
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            origin.push_back(jelly->getParticle(i));
            /*Eigen::Vector3f eulerV = jelly->getParticle(i).getAcceleration() * particleSystem.deltaTime;
            Eigen::Vector3f eulerP = jelly->getParticle(i).getVelocity() * particleSystem.deltaTime;
            jelly->getParticle(i).addVelocity(eulerV);
            jelly->getParticle(i).addPosition(eulerP); */
            jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());

        }
        particleSystem.computeJellyForce(*jelly);
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
        }
        //particleSystem.computeJellyForce(*jelly);
    }
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review ¡§ODE_basics.pptx¡¨ from p.18 - p.19

}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    for (int index = 0; index < particleSystem.jellyCount; index++) {
        auto jelly = particleSystem.getJellyPointer(index);
        std::vector<Particle> origin;
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            origin.push_back(jelly->getParticle(i));
            //Eigen::Vector3f deltaX = 0.5 * particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
            Eigen::Vector3f eulerV = jelly->getParticle(i).getAcceleration() * particleSystem.deltaTime;
            Eigen::Vector3f eulerP = jelly->getParticle(i).getVelocity() * particleSystem.deltaTime;
            jelly->getParticle(i).addVelocity(eulerV * 0.5);
            jelly->getParticle(i).addPosition(eulerP * 0.5); 
            jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());

        }
        particleSystem.computeJellyForce(*jelly);
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
            jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
        }
        //particleSystem.computeJellyForce(*jelly);
    }
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review ¡§ODE_basics.pptx¡¨ from p .18 - p .20

}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.21
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    for (int index = 0; index < particleSystem.jellyCount; index++) {
        auto jelly = particleSystem.getJellyPointer(index);
        std::vector<StateStep> k1, k2, k3, k4;
        std::vector<Particle> origin;
        for (int i = 0; i < jelly->getParticleNum(); i++)
            origin.push_back(jelly->getParticle(i));
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            StateStep k;
            k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
            k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
            jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity() / 2);
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration() / 2);
            k1.push_back(k);
            
        }
        particleSystem.computeJellyForce(*jelly);
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            StateStep k;
            k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
            k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
            jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity() / 2);
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration() / 2);
            k2.push_back(k);
        }
        particleSystem.computeJellyForce(*jelly);
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            StateStep k;
            k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
            k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
            jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
            k3.push_back(k);
        }
        particleSystem.computeJellyForce(*jelly);
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            StateStep k;
            k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
            k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
            k4.push_back(k);
        }
        for (int i = 0; i < jelly->getParticleNum(); i++) {
            jelly->getParticle(i).setPosition(origin[i].getPosition() + (k1[i].deltaPos + 2 * k2[i].deltaPos + 2 * k3[i].deltaPos + k4[i].deltaPos) / 6);
            jelly->getParticle(i).setVelocity(origin[i].getVelocity() + (k1[i].deltaVel + 2 * k2[i].deltaVel + 2 * k3[i].deltaVel + k4[i].deltaVel) / 6);
            jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
        }
    }
}
}  // namespace simulation
