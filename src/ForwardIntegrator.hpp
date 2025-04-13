/*
 * ForwardIntegrator.hpp - General forward integrator class
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicState.hpp"
#include "DynamicalSystem.hpp"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdarg>


/*
 * General forward integrator class
 */
template <class T>
class ForwardIntegrator {
public:
  ForwardIntegrator (DynamicalSystem<T>& _system, const int32_t _nstates, const int32_t _level) :
    g(_system),
    s(_nstates, DynamicState<T>(_system.getNumDerivs(), _level, 0))
  {
    // set the zero state
    s[0] = g.getInit();
  }

  // derived classes must define this method
  virtual void stepForward (const double _dt) = 0;

  // building block of many algorithms --- will only modify start if we solve for highest derivs
  DynamicState<T> EulerStep(DynamicState<T>& initial, const DynamicState<T>& derivs, const double _dt) {

    // create the state with copies of the lower derivatives
    DynamicState<T> next = initial.stepHelper();
    next.time += _dt;

    // step all lower derivatives forward
    // eg. x*  = x  + dt*x' + 0.5*dt^2*x"
    //     x'* = x' + dt*x"
    for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
      double factor = 1.0;
      for (int32_t d=nd; d>=0; d--) {
        factor *= _dt / (nd-d+1);
        next.x[d] += factor * derivs.x[nd+1];
      }
    }

    // return the state
    return next;
  }

  DynamicState<T> EulerStep(DynamicState<T>& start, const double _dt) {
    // call the general routine
    return EulerStep(start, start, _dt);
  }

  // use variadic function to merge many states into one new one
  // _in	initial state (will not use forcing terms, just state)
  // _dt	length of step to take
  // _num	number of pairs of DynamicState and coefficients passed in
  // _d1	first state of derivatives
  // _c1	first coefficient
  DynamicState<T> CombinedStep( DynamicState<T>& _in, const double _dt, const int32_t _num,
                               const DynamicState<T>& _d1, const double _c1, ...) {

    // create the new state as a copy of the initial state
    DynamicState<T> next = _in.stepHelper();
    next.time += _dt;

    //std::cout << "time: " << _in.time << "  dt: " << _dt << "\n";
    //std::cout << "state: pos " << _in.x[0].segment(0,6).transpose() << "\n";
    //std::cout << "       vel " << _in.x[1].segment(0,6).transpose() << "\n";
    //std::cout << "       acc " << _in.x[2].segment(0,6).transpose() << "\n";

    // step all lower derivatives forward using first pair _d1 and _c1
    if (_num > 0) {
      for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
        next.x[nd] += _c1*_dt * _d1.x[nd+1];
      }
    }

    // step again using additional pairs of derivatives
    va_list _derivs;
    va_start(_derivs, _c1);	// set the va pointer to the argument after _c1
    // NOTE: cannot accept references! must be full objects
    for (int32_t id=1; id<_num; ++id) {
      const DynamicState<T> deriv = va_arg(_derivs, DynamicState<T>);
      const double coef = va_arg(_derivs, double);
      // now deriv is the next DynamicState and coef is its coefficient

      // accumulate to the same state (all except top-level-derivative)
      for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
        next.x[nd] += coef*_dt * deriv.x[nd+1];
      }
    }
    va_end(_derivs);

    //std::cout << "new:   pos " << next.x[0].segment(0,6).transpose() << "\n";
    //std::cout << "       vel " << next.x[1].segment(0,6).transpose() << "\n";

    // return the state
    return next;
  }

  // use variadic function to merge many states into one new one, enhanced by higher-order derivatives!
  // for velocity systems, this is identical to CombinedStep, above!
  // _in	initial state (will not use forcing terms, just state)
  // _dt	length of step to take
  // _num	number of pairs of DynamicState and coefficients passed in
  // _d1	first state of derivatives
  // _c1	first coefficient
  DynamicState<T> CombinedStepEnhanced( DynamicState<T>& _in, const double _dt, const int32_t _num,
                                  const DynamicState<T>& _d1, const double _c1, ...) {

    // create the new state as a copy of the initial state
    DynamicState<T> next = _in.stepHelper();
    next.time += _dt;

    //std::cout << "time: " << _in.time << "  dt: " << _dt << "\n";
    //std::cout << "state: pos " << _in.x[0].segment(0,6).transpose() << "\n";
    //std::cout << "       vel " << _in.x[1].segment(0,6).transpose() << "\n";
    //std::cout << "       acc " << _in.x[2].segment(0,6).transpose() << "\n";

    // step all lower derivatives forward using first pair _d1 and _c1
    if (_num > 0) {
      for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
        double factor = _c1;
        for (int32_t d=nd; d>=0; d--) {
          factor *= _dt / (nd-d+1);
          next.x[d] += factor * _d1.x[nd+1];
        }
      }
    }

    // step again using additional pairs of derivatives
    va_list _derivs;
    va_start(_derivs, _c1);	// set the va pointer to the argument after _c1
    // NOTE: cannot accept references! must be full objects
    for (int32_t id=1; id<_num; ++id) {
      const DynamicState<T> deriv = va_arg(_derivs, DynamicState<T>);
      const double coef = va_arg(_derivs, double);
      // now deriv is the next DynamicState and coef is its coefficient

      // accumulate to the same state (all except top-level-derivative)
      for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
        double factor = coef;
        for (int32_t d=nd; d>=0; d--) {
          factor *= _dt / (nd-d+1);
          next.x[d] += factor * deriv.x[nd+1];
        }
      }
    }
    va_end(_derivs);

    //std::cout << "new:   pos " << next.x[0].segment(0,6).transpose() << "\n";
    //std::cout << "       vel " << next.x[1].segment(0,6).transpose() << "\n";

    // return the state
    return next;
  }

  T getPosition () {
    return s[0].getPos();
  }

  T getVelocity () {
    return s[0].getVel();
  }

  std::vector<T> getState () {
    return s[0].getState();
  }

  DynamicState<T> getDynamicState () {
    return s[0];
  }

  double getTime (const size_t _idx) {
    assert(_idx < s.size() && "Asking for improper index into DynamicState");
    return s[_idx].getTime();
  }

  double getTime () {
    return s[0].getTime();
  }

  T getDeriv (const int32_t deriv) {
    try {
      return s[0].x[deriv];
    } catch (std::exception& e) {
      std::cout << "Standard exception: " << e.what() << std::endl;
      return T(0);
    }
  }

  double getError (const T _trueSolution) {
    T temp = _trueSolution-getPosition();
    return( g.getErrorNorm(temp) / std::sqrt(g.getSize()) );
    //return( temp.matrix().norm() / std::sqrt(temp.size()) );
  }

protected:
  // saves the reference to the system to be integrated
  DynamicalSystem<T>& g;
  // and a series of states (current (0) and previous (1..)) to store
  std::vector<DynamicState<T>> s;
};

