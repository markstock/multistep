/*
 * DynamicalSystem.hpp - velocity and acceleration-based systems
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicState.hpp"

#include <cstdint>


/*
 * DynamicalSystem - generalized class for dynamical systems
 *
 * TODO: template this on number of derivatives?
 */
template <class T>
class DynamicalSystem {
public:
  DynamicalSystem(const int32_t _num) :
    numVars(_num)
  {}

  // number of derivatives in the dynamic state
  virtual int32_t getNumDerivs(void) = 0;
  int32_t getSize(void) { return numVars; };

  // return the initial state
  virtual DynamicState<T> getInit(void) = 0;

  virtual bool hasAccel(void) = 0;
  virtual bool hasDamping(void) = 0;
  virtual T getHighestDeriv(const T pos, const double time) = 0;
  virtual void setHighestDeriv(DynamicState<T>& state, const double time) = 0;
  virtual T getExact(const double time) = 0;
  virtual std::vector<T> getState(const double time) = 0;
  virtual double getErrorNorm(const T pos) = 0;
  virtual double getEndTime() = 0;

protected:
  // number of equations in the system
  int32_t numVars;
};


/*
 * VelocitySystem - a dynamic system driven by velocities
 *                  like vortex methods; states have x, x'
 */
template <class T>
class VelocitySystem : public DynamicalSystem<T> {
public:
  VelocitySystem(const int32_t _num) :
    DynamicalSystem<T>(_num),
    ic(DynamicState<T>(1))
  {
    ic.step = 0;
  }

  // number of derivatives in the dynamic state
  int32_t getNumDerivs(void) { return 1; }

  // return the initial state
  DynamicState<T> getInit() { return ic; }

  bool hasAccel(void) { return false; }

protected:
  // must define and store initial state
  DynamicState<T> ic;
};


/*
 * AccelerationSystem - a dynamic system driven by accelerations
 *                      like gravitation; states have x, x', x"
 */
template <class T>
class AccelerationSystem : public DynamicalSystem<T> {
public:
  AccelerationSystem(const int32_t _num) :
    DynamicalSystem<T>(_num),
    ic(DynamicState<T>(2))
  {
    ic.step = 0;
  }

  // number of derivatives in the dynamic state
  int32_t getNumDerivs(void) { return 2; }

  // return the initial state
  DynamicState<T> getInit() { return ic; }

  bool hasAccel(void) { return true; }

protected:
  // must define and store initial state
  DynamicState<T> ic;
};


/*
 * JerkSystem - a dynamic system driven by the derivative of acceleration
 *                      states have x, x', x'', x'''
 */

