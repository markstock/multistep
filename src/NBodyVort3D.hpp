/*
 * NBodyVort3D.hpp - a simple 3D vortex method
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


/*
 * A 3D vortex system with vector strengths, scalar radii
 */
class NBodyVort3D : public VelocitySystem<Eigen::ArrayXd> {
public:
  NBodyVort3D(const int32_t _num) :
    VelocitySystem<Eigen::ArrayXd>(6*_num),
    num(_num)
  {
    // num is how many bodies

    // we store the unchanging properties here (random generates numbers -1:1)
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

  // perform n-body Biot-Savart integration; uses position and strength and radius squared
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd pos, const double _time) {

    // create the accumulator vectors
    Eigen::ArrayXd newState = Eigen::ArrayXd::Zero(numVars);

    // evaluate the new vels and strengths
    // this works by performing a summation over all particles of the "influence" of the
    //   source vortex on a stick centered at the location of each target vortex, and
    //   oriented along its strength vector. this also solves for the reorientation of
    //   the strength vector
    if (true) {
      for (int32_t i=0; i<num; ++i) {
        // velocities on both "ends" of particle i
        Eigen::Vector3d thisVel(0.0, 0.0, 0.0);
        Eigen::Vector3d thisVel2(0.0, 0.0, 0.0);
        Eigen::Vector3d strVec = pos.segment(6*i+3,3).normalized();
        for (int32_t j=0; j<num; ++j) {
          // the influence of particle j
          const Eigen::Vector3d sstr = pos.segment(6*j+3,3);
          // vector between particle centers
          Eigen::Vector3d dx = pos.segment(6*j,3) - pos.segment(6*i,3);
          const double stickHalfLen = std::sqrt(radiusSquared[i]);
          // velocity at one end
          dx -= stickHalfLen*strVec;
          // eigen's .norm() returns the length of the vector, .squaredNorm() is the square of that
          double invdistsqr = 1.0 / (dx.squaredNorm()+radiusSquared(j));
          double factor = invdistsqr * std::sqrt(invdistsqr);
          thisVel[0] += factor * (dx[2]*sstr[1] - dx[1]*sstr[2]);
          thisVel[1] += factor * (dx[0]*sstr[2] - dx[2]*sstr[0]);
          thisVel[2] += factor * (dx[1]*sstr[0] - dx[0]*sstr[1]);
          // velocity at the other end
          dx += 2.0*stickHalfLen*strVec;
          invdistsqr = 1.0/(dx.norm()+radiusSquared(j));
          factor = invdistsqr * std::sqrt(invdistsqr);
          thisVel2[0] += factor * (dx[2]*sstr[1] - dx[1]*sstr[2]);
          thisVel2[1] += factor * (dx[0]*sstr[2] - dx[2]*sstr[0]);
          thisVel2[2] += factor * (dx[1]*sstr[0] - dx[0]*sstr[1]);
        }
        // these are the 3-vec velocities
        newState.segment(6*i,3) = 0.5 * (thisVel + thisVel2);
        // these are the 3-vec strengths
        newState.segment(6*i+3,3) = thisVel;
      }
    }
    return newState;
  }

  // use the best method to approximate the final state
  Eigen::ArrayXd getExact(const double _endtime) {
    int32_t maxSteps = 1000000;
    double dt = _endtime / maxSteps;
    //RK4<Eigen::ArrayXd> exact(*this,0);
    AB5<Eigen::ArrayXd> exact(*this,0,dt);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of AB5 at dt= " << dt << std::endl;
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
    // set up the return state vector, with new forces
    std::vector<Eigen::ArrayXd> state({exact.getPosition(),
                               getHighestDeriv(exact.getPosition(),_endtime)});
    return state;
  }

  // find the error norm
  double getErrorNorm(const Eigen::ArrayXd _delta) {
    return _delta.matrix().norm();
  }

  double getEndTime() {
    return 10.0;
  }

  bool hasDamping(void) { return false; }

private:
  // number of bodies
  int32_t num;
  // these values are constant in time and unique to this system
  Eigen::ArrayXd radiusSquared;
};

