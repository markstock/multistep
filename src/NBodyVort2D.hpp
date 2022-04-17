/*
 * NBodyVort2D.hpp - a simple 2D vortex method
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"
#include "MultistageIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>


/*
 * A 2D vortex system with constant strengths, radii
 */
class NBodyVort2D : public VelocitySystem<Eigen::ArrayXd> {
public:
  NBodyVort2D(const int32_t _num) :
    VelocitySystem<Eigen::ArrayXd>(2*_num),
    num(_num)
  {
    // num is how many bodies

    // we store the unchanging properties here (random generates numbers -1:1)
    circ = 2.5 * Eigen::ArrayXd::Random(num) / num;
    radiusSquared = 0.2 + 0.05 * Eigen::ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();

    // first take: distribute randomly [-1:1]
    ic.x[0] = 1.0 * Eigen::ArrayXd::Random(numVars);

    // second take: distribute in a jittered grid
    //const int32_t nx = 1 + (int32_t)std::sqrt((double)num-0.5);
    //const double dx = 2.0 / nx;
    //int32_t cnt = 0;
    //for (int32_t ix = 0; ix<nx; ++ix) {
    //for (int32_t iy = 0; iy<nx; ++iy) {
    //  if (++cnt > num) break; break;
    //}
    //}
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd pos, const double _time) {

    // create the accumulator vector
    Eigen::ArrayXd newVels = Eigen::ArrayXd::Zero(numVars);

    // evaluate the new vels
    if (true) {
      for (int32_t i=0; i<num; ++i) {
        // new velocities on particle i
        Eigen::Vector2d thisVel(0.0, 0.0);
        for (int32_t j=0; j<num; ++j) {
          // 20 flops
          // the influence of particle j
          Eigen::Vector2d dx = pos.segment(2*j,2) - pos.segment(2*i,2);
          double invdist = 1.0/(dx.norm()+radiusSquared(j));
          double factor = circ(j) * invdist * invdist;
          thisVel[0] -= dx[1] * factor;
          thisVel[1] += dx[0] * factor;
        }
        newVels.segment(2*i,2) = thisVel;
      }
    }
    return newVels;
  }

  // use the best method to approximate the final state
  Eigen::ArrayXd getExact(const double _endtime) {
    int32_t maxSteps = 100000;
    double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
    return exact.getPosition();
  }

  // return all state at the given time (here: pos and vel)
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    int32_t maxSteps = 1000;
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
  Eigen::ArrayXd circ;
  Eigen::ArrayXd radiusSquared;
};

