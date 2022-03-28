/*
 * ForwardIntegrator.hpp - General forward integrator class
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicState.hpp"
#include "DynamicalSystem.hpp"

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>


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
  DynamicState<T> EulerStep(DynamicState<T> initial, DynamicState<T> derivs, const double _dt) {

    // create the state with copies of the lower derivatives
    DynamicState<T> next = initial.stepHelper();

    // find the highest derivative - nope
    //initial.x[g.getNumDerivs()] = g.getHighestDeriv(start.x[0]);

    // step all lower derivatives forward
    for (int32_t nd=0; nd<g.getNumDerivs(); nd++) {
      double factor = 1.0;
      for (int32_t d=nd; d>=0; d--) {
        factor *= _dt / (nd-d+1);
        next.x[d] += factor * derivs.x[nd+1];
      }
    }

    // find the highest derivative - nope, again
    //next.x[g.getNumDerivs()] = g.getHighestDeriv(next.x[0]);

    // return the state
    return next;
  }

  DynamicState<T> EulerStep(DynamicState<T> start, const double _dt) {
    // call the general routine
    return EulerStep(start, start, _dt);
  }

  T getPosition () {
    return s[0].getPos();
  }

  T getVelocity () {
    return s[0].getVel();
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
    return( temp.matrix().norm() / std::sqrt(temp.size()) );
  }

protected:
  // saves the reference to the system to be integrated
  DynamicalSystem<T>& g;
  // and a series of states (current and previous) to store
  std::vector<DynamicState<T>> s;
};

