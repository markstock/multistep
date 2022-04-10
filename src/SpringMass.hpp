/*
 * SpringMass.hpp - not just about church in April
 *
 * Copyright 2022 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"

#include <cmath>

/*
 * The right way to do a sine wave, as a solution to the spring-mass system
 */
class SpringMass : public AccelerationSystem<double> {
public:
  SpringMass(const double _period) :
    AccelerationSystem<double>(1),
    period(_period)
  {
    // store initial conditions (position and velocity)
    ic.x[0] = 0.0;
    // v_0 = A sqrt(k / m) = 2 pi / period
    ic.x[1] = 2.0*M_PI/period;
  };

  // return the derivative at the given point
  double getHighestDeriv(const double _pos, const double _time) {
    // period is 2 pi sqrt(m / k)
    // acceleration is -k x / m
    // which is -x (2 pi)^2 / P^2
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    double acc = -_pos * std::pow(2.0*M_PI/period,2);
    //std::cout << "\n2nd deriv at t=" << _time << " is " << acc;
    return acc;
  }

  // just return theoretical exact position at the given time
  double getExact(const double _endtime) {
    return std::sin(2.0*M_PI*_endtime/period);
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<double> getState(const double _endtime) {
    const double pos = getExact(_endtime);
    std::vector<double> state({pos,
                               (2.0*M_PI/period) * std::cos(2.0*M_PI*_endtime/period),
                               getHighestDeriv(pos,_endtime)});
    return state;
  }

  // find the error norm
  double getErrorNorm(const double _delta) {
    //std::cout << "\n error is "
    return std::abs(_delta);
  }

protected:
  const double period;
};

