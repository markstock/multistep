/*
 * DynamicState.hpp - templated system state
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <cmath>


/*
 * A class to contain one state (one set of changing variables, and all of their derivatives)
 */
template <class T>
class DynamicState {
public:
  // delegating constructor chain
  DynamicState() :
    DynamicState(1)
  {}

  DynamicState(const int32_t _highestDeriv) :
    DynamicState(_highestDeriv, 0, 0)
  {}

  DynamicState(const int32_t _highestDeriv, const int32_t _level, const int32_t _step) :
    time(0.0),
    level(_level),
    step(_step),
    x(1+_highestDeriv, T())
  {
    // do not initialize the arrays now, wait for later
  };

  // copy constructor works just like it's supposed to

  // like a copy constructor, but advances step and zeros highest derivatives
  DynamicState stepHelper() {
    // make the new object here
    DynamicState next(x.size()-1, level, step);

    // copy value and all derivatives
    for (size_t i=0; i<x.size(); i++) {
      next.x[i] = x[i];
    }

    // highest derivative array is left as zeros
    //next.x[x.size()-1] = ArrayXd::Zero(x[0].size());
    next.x[x.size()-1] = 0.0;
    // finally, advance the step by one
    next.step++;
    // advance the "time" later
    next.time = time;

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
      std::cout << "Standard exception: " << e.what() << std::endl;
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
      std::cout << "Standard exception: " << e.what() << std::endl;
    }
  }

  T getAcc(void) {
    //if (x.size() > 2) return x[2];
    //else return ArrayXd::Zero(1);
    try {
      return x[2];
    } catch (std::exception& e) {
      std::cout << "Standard exception: " << e.what() << std::endl;
    }
  }

  // precise time of this state
  double time;
  // level 0 is at base dt, level 1 is at 2*dt, level -1 at 0.5*dt
  int32_t level;
  // step is time step within that level
  int32_t step;
  // store a value and an arbitrary number of derivatives
  // x[0] is position, x[1] velocity, x[2] acceleration, x[3] jerk, etc.
  std::vector<T> x;
};

