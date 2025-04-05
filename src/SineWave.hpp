/*
 * SineWave.hpp - just a sine wave, velocity and acceleration versions
 *
 * Copyright 2022,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"

#include <cmath>

/*
 * Two classes of sine - one where velocity is used as integration variable, one with acceleration
 */
class VelocitySine : public VelocitySystem<double> {
public:
  VelocitySine(const double _period) :
    VelocitySystem<double>(1),
    period(_period)
  {
    // store initial conditions (just position)
    ic.x[0] = 0.0;
  };

  // return the derivative at the given point
  double getHighestDeriv(const double _pos, const double _time) {
    // pos is the current value of sine, NOT time
    //std::cout << "deriv at " << _pos << " is " << std::cos(std::asin(_pos))/period << "\n";
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    double vel = std::cos(2.0*M_PI*_time/period) * (2.0*M_PI/period);
    //std::cout << "\n1st deriv at t=" << _time << " is " << vel;
    return vel;
  }

  // just return theoretical exact position at the given time
  double getExact(const double _endtime) {
    return std::sin(2.0*M_PI*_endtime/period);
  }

  // return all state at the given time (here: pos and vel)
  std::vector<double> getState(const double _endtime) {
    const double pos = getExact(_endtime);
    std::vector<double> state({pos, getHighestDeriv(pos,_endtime)});
    return state;
  }

  // find the error norm
  double getErrorNorm(const double _delta) {
    //std::cout << "\n error is "
    return std::abs(_delta);
  }

  double getEndTime() {
    return 10.0;
  }

protected:
  const double period;
};

class AccelerationSine : public AccelerationSystem<double> {
public:
  AccelerationSine(const double _period) :
    AccelerationSystem<double>(1),
    period(_period)
  {
    // store initial conditions (position and velocity)
    ic.x[0] = 0.0;
    ic.x[1] = 2.0*M_PI/period;
  };

  // return the derivative at the given point
  double getHighestDeriv(const double _pos, const double _time) {
    // pos is the current value of sine, NOT time
    //std::cout << "deriv at " << _pos << " is " << std::cos(std::asin(_pos))/period << "\n";
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    double acc = -std::sin(2.0*M_PI*_time/period) * std::pow(2.0*M_PI/period,2);
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

  double getEndTime() {
    return 10.0;
  }

protected:
  const double period;
};

