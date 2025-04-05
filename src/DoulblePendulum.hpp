/*
 * DoublePendulum.hpp - the classic chaotic system
 *
 * Copyright 2022,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"

#include <cmath>

/*
 * See https://scienceworld.wolfram.com/physics/DoublePendulum.html
 *
 * Inputs are masses and initial positions, vels are assumed zero
 * Two primary variables: angular positions
 * NOT DONE
 */
class DoublePendulum : public AccelerationSystem<std::array<double,2>> {
public:
  DoublePendulum(const double _th1, const double _m1, const double _th2, const double _m2) :
    AccelerationSystem<double>(2),
    m1(_m1), m2(_m2), th1(_th1), th2(_th2)
  {
    // store initial conditions (position and velocity)
    ic.x[0] = 0.0;
    // v_0 = A sqrt(k / m) = 2 pi / period
    ic.x[1] = 2.0*M_PI/period;
  };

  // return the derivative at the given point
  std::array<double,2> getHighestDeriv(const std::array<double,2> _pos, const double _time) {
    // period is 2 pi sqrt(m / k)
    // acceleration is -k x / m
    // which is -x (2 pi)^2 / P^2
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    std::array<double,2> acc = -_pos * std::pow(2.0*M_PI/period,2);
    //std::cout << "\n2nd deriv at t=" << _time << " is " << acc;
    return acc;
  }

  // just return theoretical exact position at the given time
  std::array<double,2> getExact(const double _endtime) {
    return std::sin(2.0*M_PI*_endtime/period);
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<std::array<double,2>> getState(const double _endtime) {
    const std::array<double,2> pos = getExact(_endtime);
    std::vector<std::array<double,2>> state({pos,
                               (2.0*M_PI/period) * std::cos(2.0*M_PI*_endtime/period),
                               getHighestDeriv(pos,_endtime)});
    return state;
  }

  // find the error norm
  double getErrorNorm(const std::array<double,2> _delta) {
    //std::cout << "\n error is "
    return std::abs(_delta);
  }

  double getEndTime() {
    return 10.0;
  }

protected:
  const double m1, m2, th1, th2;
};

