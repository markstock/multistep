/*
 * LennardJones.hpp - anharmonic oscillator
 *
 * Copyright 2022 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"

#include <cmath>

/*
 * The right way to do a sine wave, as a solution to the spring-mass system
 */
class LennardJones : public AccelerationSystem<double> {
public:
  LennardJones(const double _xs, const double _mass) :
    AccelerationSystem<double>(1),
    xstart(_xs), mass(_mass)
  {
    // store initial conditions (position and velocity)
    ic.x[0] = xstart;
    ic.x[1] = 0.0;
  };

  // return the derivative at the given point
  double getHighestDeriv(const double _pos, const double _time) {
    double acc = (std::pow(_pos,-12) - std::pow(_pos,-6)) / mass;
    //std::cout << "\n2nd deriv at x=" << _pos << " and t=" << _time << " is " << acc;
    return acc;
  }

  // just return theoretical exact position at the given time
  double getExact(const double _endtime) {
    int32_t maxSteps = 1000000;
    double dt = _endtime / maxSteps;
    RK4<double> exact(*this,0);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) { exact.stepForward(dt); }
    return exact.getPosition();
  }

  // return all state at the given time (here: pos, vel, acc),
  // this is to populate back history for the multi-step integrators
  std::vector<double> getState(const double _endtime) {
    int32_t maxSteps = 10000;
    double dt = _endtime / maxSteps;
    RK4<double> exact(*this,0);
    //std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) { exact.stepForward(dt); }

    // set up the return state vector, with new forces
    std::vector<double> state({exact.getPosition(),
                               exact.getVelocity(),
                               getHighestDeriv(exact.getPosition(),_endtime)});
    return state;
  }

  // find the error norm
  double getErrorNorm(const double _delta) {
    //std::cout << "\n error is "
    return std::abs(_delta);
  }

protected:
  const double xstart, mass;
};

