/*
 * NBodyGrav3D.hpp - a simple 3D gravitational system
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
 * A gravitational n-body system
 */
class NBodyGrav3D : public AccelerationSystem<Eigen::ArrayXd> {
public:
  NBodyGrav3D(const int32_t _num) :
    AccelerationSystem<Eigen::ArrayXd>(3*_num),
    num(_num)
  {
    // num is how many bodies
    if (_num == 5) {
      // we are doing the comuter language benchmark simulation of the 5 massive bodies in the solar system
      // ic.x[0] is vector of initial positions: 3 values per body
      // ic.x[1] is vector of initial velocities: 3 values per body
      ic.x[0] = Eigen::ArrayXd::Zero(3*num);
      ic.x[1] = Eigen::ArrayXd::Zero(3*num);
      mass = Eigen::ArrayXd::Zero(num);

      // sun
      ic.x[0].segment(0,3) << 0.0, 0.0, 0.0;
      ic.x[1].segment(0,3) << 0.0, 0.0, 0.0;
      mass[0] = 1.0;
      size_t ib = 0;

      // jupiter
      ib += 3;
      ic.x[0].segment(3,3) << 4.84143144246472090e+00, -1.16032004402742839e+00, -1.03622044471123109e-01;
      ic.x[1].segment(3,3) << 1.66007664274403694e-03, 7.69901118419740425e-03, -6.90460016972063023e-05;
      mass[1] = 9.54791938424326609e-04;

      // saturn
      ib += 3;
      ic.x[0].segment(6,3) << 8.34336671824457987e+00, 4.12479856412430479e+00, -4.03523417114321381e-01;
      ic.x[1].segment(6,3) << -2.76742510726862411e-03, 4.99852801234917238e-03, 2.30417297573763929e-05;
      mass[2] = 2.85885980666130812e-04;

      // uranus
      ib += 3;
      ic.x[0].segment(9,3) << 1.28943695621391310e+01, -1.51111514016986312e+01, -2.23307578892655734e-01;
      ic.x[1].segment(9,3) << 2.96460137564761618e-03, 2.37847173959480950e-03, -2.96589568540237556e-05;
      mass[3] = 4.36624404335156298e-05;

      // neptune
      ib += 3;
      ic.x[0].segment(12,3) << 1.53796971148509165e+01, -2.59193146099879641e+01, 1.79258772950371181e-01;
      ic.x[1].segment(12,3) << 2.68067772490389322e-03, 1.62824170038242295e-03, -9.51592254519715870e-05;
      mass[4] = 5.15138902046611451e-05;

      const double pi = 3.141592653589793;
      const double days_per_year = 365.24;
      ic.x[1] = ic.x[1] * days_per_year;
      mass = mass * 4.0*pi*pi;

      // exact integration: no smoothing radius
      radiusSquared = Eigen::ArrayXd::Zero(num);

    } else if (_num == 7) {
      // we are doing the Pleiades problem
      // ic.x[0] is vector of initial positions: 3 values per body
      // ic.x[1] is vector of initial velocities: 3 values per body
      ic.x[0] = Eigen::ArrayXd::Zero(3*num);
      ic.x[1] = Eigen::ArrayXd::Zero(3*num);
      mass = Eigen::ArrayXd::Zero(num);

      // the bodies
      ic.x[0].segment(0,3) << 3.0, 3.0, 0.0;
      ic.x[1].segment(0,3) << 0.0, 0.0, 0.0;
      mass[0] = 1.0;

      ic.x[0].segment(3,3) << 3.0, -3.0, 0.0;
      ic.x[1].segment(3,3) << 0.0, 0.0, 0.0;
      mass[1] = 2.0;

      ic.x[0].segment(6,3) << -1.0, 2.0, 0.0;
      ic.x[1].segment(6,3) << 0.0, 0.0, 0.0;
      mass[2] = 3.0;

      ic.x[0].segment(9,3) << -3.0, 0.0, 0.0;
      ic.x[1].segment(9,3) << 0.0, -1.25, 0.0;
      mass[3] = 4.0;

      ic.x[0].segment(12,3) << 2.0, 0.0, 0.0;
      ic.x[1].segment(12,3) << 0.0, 1.0, 0.0;
      mass[4] = 5.0;

      ic.x[0].segment(15,3) << -2.0, -4.0, 0.0;
      ic.x[1].segment(15,3) << 1.75, 0.0, 0.0;
      mass[5] = 6.0;

      ic.x[0].segment(18,3) << 2.0, 4.0, 0.0;
      ic.x[1].segment(18,3) << -1.5, 0.0, 0.0;
      mass[6] = 7.0;

      // exact integration: no smoothing radius
      radiusSquared = Eigen::ArrayXd::Zero(num);

    } else {
      // just generate random bodies instead

      // we store the unchanging properties here
      mass = (2.0 + Eigen::ArrayXd::Random(num)) / (2.0*num);
      radiusSquared = 0.2 + 0.1 * Eigen::ArrayXd::Random(num);
      radiusSquared = radiusSquared.square();

      // store initial conditions in case we want to reuse this
      ic.x[0] = 10.0 * Eigen::ArrayXd::Random(numVars);
      ic.x[1] = 1.0 * Eigen::ArrayXd::Random(numVars);
      //std::cout << "NBodyGrav3D::NBodyGrav3D " << ic.x[0].segment(0,4).transpose() << std::endl;
    }

    // offset the momentum by assigning a velocity to the 0th mass
    Eigen::ArrayXd momentum = Eigen::ArrayXd::Zero(3);
    for (int32_t i=1; i<num; ++i) {
      momentum += mass[i] * ic.x[1].segment(3*i,3);
    }
    ic.x[1].segment(0,3) -= momentum / mass[0];
    //for (int32_t i=1; i<3; ++i) ic.x[1][i] = momentum[i] / mass[0];
    if (num == 5) std::cout << "NBodyGrav3D set sun vel to " << ic.x[1].segment(0,3).transpose() << std::endl;

    // and calculate initial energy
    initial_energy = getEnergy(ic);
    if (num == 5) std::cout << "  starting energy " << std::setprecision(16) << initial_energy << std::endl;
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd pos, const double _time) {

    // generate the output vector (3 accelerations per projectile)
    Eigen::ArrayXd newVal = Eigen::ArrayXd::Zero(numVars);

    // evaluate forces
    for (int32_t i=0; i<num; ++i) {
      // new accelerations on particle i
      Eigen::Vector3d newAcc(0.0, 0.0, 0.0);
      for (int32_t j=0; j<num; ++j) {
        if (i != j) {
          // 20 flops
          // the influence of particle j
          Eigen::Vector3d dx = pos.segment(3*j,3) - pos.segment(3*i,3);
          double invdist = 1.0/(dx.norm()+radiusSquared(j));
          newAcc += dx * (mass(j) * invdist * invdist * invdist);
        }
      }
      newVal.segment(3*i,3) = newAcc;
    }
    return newVal;
  }

  // use the best method to approximate the final state
  // initial is -0.1690751638285245
  //   RK4 gets -0.1690751638285206
  //  verlet is -0.1471901297585013 ?!?
  //   AB5 gets -0.1690751638285266 - closer!
  Eigen::ArrayXd getExact(const double _endtime) {
    int32_t maxSteps = 10000000;
    double dt = _endtime / maxSteps;
    //RK4<Eigen::ArrayXd> exact(*this,0);
    //Verlet<Eigen::ArrayXd> exact(*this,0,dt);
    AB5<Eigen::ArrayXd> exact(*this,0,dt);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of AB5 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
      //if (i%10000 == 0) std::cout << "    " << dt*i << "  " << exact.getPosition().segment(3,3).transpose() << std::endl;
    }
    if (num == 5) std::cout << "  ending energy " << std::setprecision(16) << getEnergy(exact.getDynamicState()) << std::endl;
    return exact.getPosition();
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    int32_t maxSteps = 10000;
    double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
    //return exact.getState();
    // set up the return state vector, with new forces
    std::vector<Eigen::ArrayXd> state({exact.getPosition(),
                               exact.getVelocity(),
                               getHighestDeriv(exact.getPosition(),_endtime)});
    return state;
  }

  // use state to calculate energy
  double getEnergy(const DynamicState<Eigen::ArrayXd>& _state) {
    double energy = 0.0;

    // add kinetic energy
    for (int32_t i=0; i<num; ++i) {
      const Eigen::Vector3d vel = _state.x[1].segment(3*i,3);
      energy += 0.5 * mass[i] * vel.dot(vel);
    }

    // add gravitational potential energy from every pair
    for (int32_t i=0; i<num; ++i) {
      for (int32_t j=i+1; j<num; ++j) {
        Eigen::Vector3d dx = _state.x[0].segment(3*j,3) - _state.x[0].segment(3*i,3);
        double invdist = 1.0/dx.norm();
        energy -= mass[i] * mass[j] * invdist;
      }
    }

    return energy;
  }

  // find the error norm
  double getErrorNorm(const Eigen::ArrayXd _delta) {
    return _delta.matrix().norm();
  }

private:
  // number of bodies
  int32_t num;
  // these values are constant in time and unique to this system
  Eigen::ArrayXd mass;
  Eigen::ArrayXd radiusSquared;
  // save energy to compare
  double initial_energy = 0.0;
};

