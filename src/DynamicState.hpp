/*
 * DynamicState.hpp - templated system state
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>
#include <span>

// for potential use in output, derivName[0] being position (not a derivative, really, hence 0)
std::vector<std::string> derivName{ "position", "velocity", "acceleration", "jerk", "snap", "crackle", "pop" };

/*
 * A class to contain one state (one set of changing variables, and all of their derivatives)
 */
template <class T>
class DynamicState {
public:
  // delegating constructor chain
  DynamicState() :
    // defaults to system with x and x' only
    DynamicState(1)
  {}

  DynamicState(const int32_t _highestDeriv) :
    DynamicState(_highestDeriv, 0, 0)
  {}

  DynamicState(const int32_t _highestDeriv, const int32_t _level, const int32_t _step) :
    DynamicState(_highestDeriv, 0.0, _level, _step)
  {}

  DynamicState(const int32_t _highestDeriv, const double _time, const int32_t _level, const int32_t _step) :
    time(_time),
    level(_level),
    step(_step),
    x(1+_highestDeriv, T()),
    tldCurrent(false)
  {
    // do not initialize the arrays now, wait for later
  };

  // copy constructor works just like it's supposed to

  // like a copy constructor, but advances step and zeros highest derivatives
  DynamicState stepHelper() {
    // make the new object here
    // advance the step by one, but keep the time
    DynamicState next(x.size()-1, time, level, step+1);

    // copy value and all derivatives
    for (size_t i=0; i<x.size(); i++) {
      next.x[i] = x[i];
    }

    // highest derivative array is left as zeros
    //next.x[x.size()-1] = ArrayXd::Zero(x[0].size());
    next.x[x.size()-1] = 0.0;
    tldCurrent = false;

    return next;
  }

  double getTimeStep () {
    return step * std::pow(2.0, level);
  }

  double getTime () {
    return time;
  }

  T getPos(void) {
    //std::cout << "DynamicState::getPos " << x.size() << endl;
    //std::cout << "DynamicState::getPos " << x.size() << "  " << x[0].size() << endl;
    //std::cout << "DynamicState::getPos " << x.size() << "  " << x[0].size() << "  " << x[0].segment(0,4).transpose() << endl;
    try {
      return x[0];
    } catch (std::exception& e) {
      std::cout << "No position: " << e.what() << std::endl;
      return T(0);
    }
    //if (x.size() > 0) return x[0];
    //else return ArrayXd::Zero(1);
  }

  T getVel(void) {
    //if (x.size() > 1) return x[1];
    //else return ArrayXd::Zero(1);
    try {
      return x[1];
    } catch (std::exception& e) {
      std::cout << "No velocity: " << e.what() << std::endl;
      exit(1);
    }
  }

  T getAcc(void) {
    //if (x.size() > 2) return x[2];
    //else return ArrayXd::Zero(1);
    try {
      return x[2];
    } catch (std::exception& e) {
      std::cout << "No acceleration: " << e.what() << std::endl;
      exit(1);
    }
  }

  T getDeriv(const int32_t _deriv) {
    try {
      return x[_deriv];
    } catch (std::exception& e) {
      std::cout << "No " << _deriv << " derivative: " << e.what() << std::endl;
      exit(1);
    }
  }

  // old: returns the entire vector
  std::vector<T> getState(void) {
    return x;
  }

  // new:
  // state is every value up to the top-level derivative (x and x' for acceleration systems)
  std::span<T> getStateVec(void) {
    return std::span<T>{x.data(), x.size()-1};
  }
  // deriv is every value above the positions (x' and x'' for acceleration systems)
  std::span<T> getDerivVec(void) {
    return std::span<T>{x.data()+1, x.size()-1};
  }

  // precise time of this state
  double time;
  // level 0 is at base dt, level 1 is at 2*dt, level -1 at 0.5*dt
  // this is only useful for states that are part of an adaptive time-stepping scheme
  int32_t level;
  // step is time step at level 0
  int32_t step;
  // store a value and an arbitrary number of derivatives
  // x[0] is position, x[1] velocity, x[2] acceleration, x[3] jerk, etc.
  std::vector<T> x;
  // need to know if top-level derivatives are current with the rest of the state
  bool tldCurrent;
};

