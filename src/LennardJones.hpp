/*
 * LennardJones.hpp - anharmonic oscillator
 *
 * Copyright 2022,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"
#include "MultistageIntegrator.hpp"
#include "MultistepIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <cmath>

/*
 * A system of a number of particles (electrons) under Lennard-Jones forces (potential?)
 */
class LennardJones : public AccelerationSystem<Eigen::ArrayXd> {
public:
  LennardJones(const int32_t _num, const double _xs, const double _mass) :
    AccelerationSystem<Eigen::ArrayXd>(_num),
    num(_num), xstart(_xs), avgmass(_mass)
  {
    // set unique masses
    mass = avgmass * (1.0 + 0.1*Eigen::ArrayXd::Random(num));

    // store initial conditions (position and velocity)
    ic.x[0] = 1.5 + 0.5*Eigen::ArrayXd::Random(num);
    ic.x[1] = Eigen::ArrayXd::Zero(num);

    //std::cout << "initial particles are:\n";
    //for (int32_t i=0; i<num; ++i) std::cout << "\t" << i << "\t" << mass[i] << "\t" << ic.x[0][i] << "\n";
  };

  // return the derivative at the given point
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd _pos, const double _time) {
    Eigen::ArrayXd acc = (_pos.pow(-12) - _pos.pow(-6)) / mass;
    //std::cout << "\n2nd deriv at x=" << _pos << " and t=" << _time << " is " << acc;
    return acc;
  }

  // set the derivative at the given point
  void setHighestDeriv(DynamicState<Eigen::ArrayXd>& _state, const double _time) {
    const Eigen::ArrayXd& pos = _state.x[0];
    Eigen::ArrayXd& acc = _state.x[2];
    acc = (pos.pow(-12) - pos.pow(-6)) / mass;
    //std::cout << "\n2nd deriv at x=" << _pos << " and t=" << _time << " is " << acc;
    return;
  }

  // just return theoretical exact position at the given time
  Eigen::ArrayXd getExact(const double _endtime) {
    const int32_t maxSteps = 1000000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) { exact.stepForward(dt); }
    return exact.getPosition();
  }

  // return all state at the given time (here: pos, vel, acc),
  // this is to populate back history for the multi-step integrators
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    const int32_t maxSteps = 10000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    for (int32_t i=0; i<maxSteps; ++i) { exact.stepForward(dt); }

    // set up the return state vector, with new forces
    std::vector<Eigen::ArrayXd> state({exact.getPosition(),
                               exact.getVelocity(),
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

protected:
  // number of bodies
  int32_t num;
  // these values are the averages for the auto-generated points
  const double xstart, avgmass;
  Eigen::ArrayXd mass;
};

