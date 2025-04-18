/*
 * SpringMass.hpp - not just about church in April
 *
 * Copyright 2022,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicalSystem.hpp"
#include "MultistageIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <cmath>

/*
 * The right way to do a sine wave, as a solution to the spring-mass system
 */
class SpringMass : public AccelerationSystem<Eigen::ArrayXd> {
public:
  SpringMass(const int32_t _num, const double _period) :
    AccelerationSystem<Eigen::ArrayXd>(_num),
    num(_num), avgperiod(_period)
  {
    // set unique periods
    period = avgperiod * (1.0 + 0.1*Eigen::ArrayXd::Random(num));

    // store initial conditions (position and velocity)
    ic.x[0] = Eigen::ArrayXd::Zero(num);
    // v_0 = A sqrt(k / m) = 2 pi / period
    ic.x[1] = 2.0*M_PI/period;
  };

  // return the derivative at the given point
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd _pos, const double _time) {
    // period is 2 pi sqrt(m / k)
    // acceleration is -k x / m
    // which is -x (2 pi)^2 / P^2
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    const Eigen::ArrayXd tmp = 2.0*M_PI / period;
    Eigen::ArrayXd acc = -_pos * tmp.pow(2);
    //std::cout << "\n2nd deriv at t=" << _time << " is " << acc;
    return acc;
  }

  // set the derivative at the given point
  void setHighestDeriv(DynamicState<Eigen::ArrayXd>& _state) {
    const Eigen::ArrayXd& pos = _state.x[0];
    Eigen::ArrayXd& acc = _state.x[2];
    const Eigen::ArrayXd tmp = 2.0*M_PI / period;
    acc = -pos * tmp.pow(2);
    return;
  }

  // just return theoretical exact position at the given time
  Eigen::ArrayXd getExact(const double _endtime) {
    return Eigen::sin(2.0*M_PI*_endtime/period);
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    const Eigen::ArrayXd pos = getExact(_endtime);
    const Eigen::ArrayXd tmp = 2.0*M_PI / period;
    // no need to integrate with fine resolution, we know the answer already
    std::vector<Eigen::ArrayXd> state({pos,
                               tmp * Eigen::cos(tmp*_endtime),
                               getHighestDeriv(pos,_endtime)});
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
  const double avgperiod;
  Eigen::ArrayXd period;
  //Eigen::ArrayXd xstart;
};


/*
 * SpringMassDamper - add a velocity term to the 1-D equations
 */
class SpringMassDamper : public AccelerationSystem<Eigen::ArrayXd> {
public:
  SpringMassDamper(const int32_t _num, const double _period) :
    AccelerationSystem<Eigen::ArrayXd>(_num),
    num(_num), avgperiod(_period)
  {
    // set unique periods
    period = avgperiod * (1.0 + 0.1*Eigen::ArrayXd::Random(num));

    // and varying damping
    damping = 0.1 + 0.05*Eigen::ArrayXd::Random(num);

    // store initial conditions (position and velocity)
    ic.x[0] = Eigen::ArrayXd::Zero(num);
    // v_0 = A sqrt(k / m) = 2 pi / period
    ic.x[1] = 2.0*M_PI/period;
  };

  // return the derivative at the given point
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd _pos, const double _time) {
    assert(false and "SpringMassDamper::getHighestDeriv impossible!");
    // period is 2 pi sqrt(m / k)
    // acceleration is -k x / m
    // which is -x (2 pi)^2 / P^2
    //return std::cos(std::asin(_pos))/period;
    // no, must use current time to determine this
    const Eigen::ArrayXd tmp = 2.0*M_PI / period;
    Eigen::ArrayXd acc = -_pos * tmp.pow(2);
    //Eigen::ArrayXd acc = -_pos * tmp.pow(2) - damping * ;
    //std::cout << "\n2nd deriv at t=" << _time << " is " << acc;
    return acc;
  }

  // set the derivative at the given point
  void setHighestDeriv(DynamicState<Eigen::ArrayXd>& _state) {
    const Eigen::ArrayXd& pos = _state.x[0];
    const Eigen::ArrayXd& vel = _state.x[1];
    Eigen::ArrayXd& acc = _state.x[2];
    const Eigen::ArrayXd tmp = 2.0*M_PI / period;
    acc = -pos * tmp.pow(2) - vel*damping;
    return;
  }

  // no theoretical exact position, so simulate it
  Eigen::ArrayXd getExact(const double _endtime) {
    // run getState(_endtime) and return just the velocity portion
    //const std::vector<Eigen::ArrayXd> endstate = getState(_endtime);
    //return endstate[0];

    const int32_t maxSteps = 1000000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
    //std::cout << "Exact positions:\n" << exact.getPosition() << std::endl;
    return exact.getPosition();
  }

  // return all state at the given time (here: pos, vel, acc)
  std::vector<Eigen::ArrayXd> getState(const double _endtime) {
    const int32_t maxSteps = 10000;
    const double dt = _endtime / maxSteps;
    RK4<Eigen::ArrayXd> exact(*this,0);
    //std::cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << std::endl;
    for (int32_t i=0; i<maxSteps; ++i) {
      exact.stepForward(dt);
    }
    // set up the return state vector, with new forces
    DynamicState<Eigen::ArrayXd> state = exact.getDynamicState();
    setHighestDeriv(state);
    return state.x;
  }

  // find the error norm
  double getErrorNorm(const Eigen::ArrayXd _delta) {
    return _delta.matrix().norm();
  }

  double getEndTime() {
    return 10.0;
  }

  bool hasDamping(void) { return true; }

protected:
  // number of bodies
  int32_t num;
  const double avgperiod;
  Eigen::ArrayXd period;
  Eigen::ArrayXd damping;
  //Eigen::ArrayXd xstart;
};

