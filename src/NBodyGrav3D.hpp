/*
 * NBodyGrav3D.hpp - a simple 3D gravitational system
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#include "DynamicalSystem.hpp"
#include "MultistageIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>


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

    // we store the unchanging properties here
    mass = (2.0 + Eigen::ArrayXd::Random(num)) / (2.0*num);
    radiusSquared = 0.2 + 0.1 * Eigen::ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();

    // store initial conditions in case we want to reuse this
    ic.x[0] = 10.0 * Eigen::ArrayXd::Random(numVars);
    ic.x[1] = 1.0 * Eigen::ArrayXd::Random(numVars);
    //std::cout << "NBodyGrav3D::NBodyGrav3D " << ic.x[0].segment(0,4).transpose() << std::endl;
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd pos, const double _time) {

    // generate the output vector
    Eigen::ArrayXd newVal = Eigen::ArrayXd::Zero(numVars);

    // evaluate forces
    if (true) {
      for (int32_t i=0; i<num; ++i) {
        // new accelerations on particle i
        Eigen::Vector3d newAcc(0.0, 0.0, 0.0);
        for (int32_t j=0; j<num; ++j) {
          // 20 flops
          // the influence of particle j
          Eigen::Vector3d dx = pos.segment(3*j,3) - pos.segment(3*i,3);
          double invdist = 1.0/(dx.norm()+radiusSquared(j));
          newAcc += dx * (mass(j) * invdist * invdist * invdist);
        }
        newVal.segment(3*i,3) = newAcc;
      }
    }
    return newVal;
  }

  // use the best method to approximate the final state
  Eigen::ArrayXd getExact(const double _endtime) {
    int32_t maxSteps = 10000;
    double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
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
    return exact.getState();
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
};

