/*
 * Projectiles3D.hpp - many projectiles in 3D with drag
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"
#include "MultistageIntegrator.hpp"
#include "MultistepIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>
#include <iomanip>


/*
 * A set of independent projectiles
 */
class Projectiles3D : public AccelerationSystem<Eigen::ArrayXd> {
public:
  Projectiles3D(const int32_t _num) :
    AccelerationSystem<Eigen::ArrayXd>(3*_num),
    num(_num)
  {
    // num is how many bodies
    // ic.x[0] is vector of initial positions: 3 values per body
    // ic.x[1] is vector of initial velocities: 3 values per body

    // initialize the arrays
    ic.x[0] = Eigen::ArrayXd::Zero(3*num);
    ic.x[1] = Eigen::ArrayXd::Zero(3*num);
    mass = Eigen::ArrayXd::Zero(num);
    area = Eigen::ArrayXd::Zero(num);
    cD = Eigen::ArrayXd::Zero(num);

    // we store the unchanging properties here
    // mass is initially something like density (kg / dl)
    mass = 2.3 + 2. * Eigen::ArrayXd::Random(num);
    // area is from bullet to basketball
    area = 0.01 + 0.3 * Eigen::ArrayXd::Random(num);
    area = area.square();
    // cD is from bullet to brick
    cD = Eigen::ArrayXd::Random(num);
    cD = 0.005 + cD.square();
    // now adjust the mass to something more reasonable
    mass = mass * area * 100. * cD;

    // projectile elevation angle (theta) and orientation (phi) and initial speed
    Eigen::ArrayXd theta = 0.5 + 0.4 * Eigen::ArrayXd::Random(num);
    Eigen::ArrayXd phi = 3.1416 * Eigen::ArrayXd::Random(num);
    Eigen::ArrayXd vel = 30. + 20. * Eigen::ArrayXd::Random(num);
    //ic.x[1] = 1.0 * Eigen::ArrayXd::Random(numVars);
    // set the vector velocities
    for (int32_t i=0; i<num; ++i) {
      ic.x[1][3*i+0] = vel[i] * std::cos(theta[i]) * std::cos(phi[i]);
      ic.x[1][3*i+1] = vel[i] * std::cos(theta[i]) * std::sin(phi[i]);
      ic.x[1][3*i+2] = vel[i] * std::sin(theta[i]);
      //std::cout << "projectile mass " << mass[i] << " area " << area[i] << " cd " << cD[i] << " and vel " << ic.x[1].segment(3*i,3).transpose() << std::endl;
    }
    //std::cout << "Projectiles3D::Projectiles3D " << ic.x[0].segment(0,4).transpose() << std::endl;
  };

  // old way: cannot compute forces without velocities
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd _pos, const double _time) {
    assert(false and "getHighestDeriv incomplete!");
    return Eigen::ArrayXd::Zero(numVars);
  }

  // perform force calculation; uses velocity, mass, cD (but not position!)
  void setHighestDeriv(DynamicState<Eigen::ArrayXd>& _state, const double _time) {
    //const Eigen::ArrayXd& pos = _state.x[0];
    const Eigen::ArrayXd& vel = _state.x[1];
    Eigen::ArrayXd& acc = _state.x[2];

    // generate the output vector (3 accelerations per projectile)
    acc = Eigen::ArrayXd::Zero(numVars);

    // evaluate forces
    for (int32_t i=0; i<num; ++i) {
      // new accelerations on particle i
      Eigen::Vector3d newAcc(0.0, 0.0, grav);
      Eigen::Vector3d thisVel = vel.segment(3*i,3);
      const double drag = 0.5 * cD[i] * area[i] * rho * thisVel.squaredNorm();
      // normalize velocity to find direction of drag force
      thisVel.norm();
      newAcc -= thisVel * (drag / mass[i]);
      acc.segment(3*i,3) = newAcc;
      //std::cout << "time " << _time << " projectile " << i << " has vel " << vel.segment(3*i,3).transpose() << "\n";
    }
    return;
  }

  // use the best method to approximate the final state
  Eigen::ArrayXd getExact(const double _endtime) {
    const int32_t maxSteps = 1000000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    //Verlet<Eigen::ArrayXd> exact(*this,0,dt);
    //AB5<Eigen::ArrayXd> exact(*this,0,dt);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
      //if (i%2000 == 0) std::cout << "    " << dt*i << "  " << exact.getPosition().segment(3,3).transpose() << std::endl;
    }
    //std::cout << "Exact positions:\n" << exact.getPosition() << std::endl;
    return exact.getPosition();
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    const int32_t maxSteps = 10000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
    // set up the return state vector, with new forces
    DynamicState<Eigen::ArrayXd> state = exact.getDynamicState();
    setHighestDeriv(state, _endtime);
    return state.x;
  }

  // find the error norm
  double getErrorNorm(const Eigen::ArrayXd _delta) {
    return _delta.matrix().norm();
  }

  double getEndTime() {
    return 4.0;
  }

  bool hasDamping(void) { return true; }

private:
  // number of bodies
  int32_t num;
  // these values are constant in time and unique to this system
  Eigen::ArrayXd mass;				// kg
  Eigen::ArrayXd area;				// in m^2
  Eigen::ArrayXd cD;				// non-dimensional
  // physical costants
  const double grav = -9.80665;		// m / s^2
  const double rho = 1.2754;		// kg / m^3
};

