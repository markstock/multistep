/*
 * MultistepIntegrator.hpp - Multi-step forward integrator classes (using history)
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicState.hpp"
#include "DynamicalSystem.hpp"
#include "ForwardIntegrator.hpp"
#include "MultistageIntegrator.hpp"

#include <cstdint>
#include <cassert>
#include <iostream>


/*
 * Multi-step forward integrator class (saves previous solutions)
 */
template <class T>
class MultistepIntegrator : public ForwardIntegrator<T> {
public:
  MultistepIntegrator (const int32_t _nsteps, DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    ForwardIntegrator<T>(_system, _nsteps, _level)
  {
    // the t=0 state is set in parent constructor
    // previous states are set here by calling the system's getExact function
    for (int32_t istep=1; istep<_nsteps; ++istep) {
      const double thistime = -1.0 * _dt * istep;
      this->s[istep].step = -istep;
      this->s[istep].time = thistime;
      this->s[istep].x = this->g.getState(thistime);
    }
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 2nd order
 */
template <class T>
class AB2 : public MultistepIntegrator<T> {
public:
  AB2 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(2, _system, _level, _dt)
  {}
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    // assert that we have at least two states in the history
    //static_assert(s.size() >= 2, "State vector does not have enough entries");
    // assert that this dt matches the previous dt
    const int32_t numDerivs = this->g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    //this->s[0].x[numDerivs] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration: AB2 for first one down, AM2 for all others
    this->s[0].x[numDerivs-1] = this->s[1].x[numDerivs-1] + 0.5 * _dt * (3.0*this->s[1].x[numDerivs] - this->s[2].x[numDerivs]);
    for (int32_t deriv=numDerivs-1; deriv>0; deriv--) {
      this->s[0].x[deriv-1] += 0.5 * _dt * (this->s[0].x[deriv] + this->s[1].x[deriv]);
    }

    // perform forward integration: AB2 for velocity, AM2 for position
    //this->s[0].x[1] = this->s[1].x[1] + 0.5 * _dt * (3.0*this->s[1].x[2] - this->s[2].x[2]);
    //this->s[0].x[0] = this->s[1].x[0] + 0.5 * _dt * (this->s[0].x[1] + this->s[1].x[1]);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 3rd order
 */
template <class T>
class AB3 : public MultistepIntegrator<T> {
public:
  AB3 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(3, _system, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    const int32_t nd = this->g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration: AB3 for first, AM3 for all lower-order derivatives
    this->s[0].x[nd-1] += _dt * (23.0*this->s[1].x[nd] - 16.0*this->s[2].x[nd] + 5.0*this->s[3].x[nd]) / 12.0;
    for (int32_t deriv=nd-1; deriv>0; deriv--) {
      this->s[0].x[deriv-1] += _dt * (5.0*this->s[0].x[deriv] + 8.0*this->s[1].x[deriv] - this->s[2].x[deriv]) / 12.0;
    }

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 4th order
 */
template <class T>
class AB4 : public MultistepIntegrator<T> {
public:
  AB4 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(4, _system, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    const int32_t nd = this->g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration: AB4 for first, AM4 for all lower-order derivatives
    this->s[0].x[nd-1] += _dt * (55.0*this->s[1].x[nd] - 59.0*this->s[2].x[nd] + 37.0*this->s[3].x[nd] - 9.0*this->s[4].x[nd]) / 24.0;
    for (int32_t deriv=nd-1; deriv>0; deriv--) {
      this->s[0].x[deriv-1] += _dt * (9.0*this->s[0].x[deriv] + 19.0*this->s[1].x[deriv] - 5.0*this->s[2].x[deriv] + this->s[3].x[deriv]) / 24.0;
    }

    // perform forward integration: AB4 for velocity, AM4 for position
    //this->s[0].x[1] = this->s[1].x[1] +  _dt * (55.0*this->s[1].x[2] - 59.0*this->s[2].x[2] + 37.0*this->s[3].x[2] - 9.0*this->s[4].x[2]) / 24.0;
    //this->s[0].x[0] = this->s[1].x[0] +  _dt * (9.0*this->s[0].x[1] + 19.0*this->s[1].x[1] - 5.0*this->s[2].x[1] + this->s[3].x[1]) / 24.0;

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 5th order
 */
template <class T>
class AB5 : public MultistepIntegrator<T> {
public:
  AB5 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(5, _system, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    const int32_t nd = this->g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration: AB5 for first, AM5 for all lower-order derivatives
    this->s[0].x[nd-1] += _dt * (1901.0*this->s[1].x[nd] - 2774.0*this->s[2].x[nd] + 2616.0*this->s[3].x[nd] - 1274.0*this->s[4].x[nd] + 251.0*this->s[5].x[nd]) / 720.0;
    for (int32_t deriv=nd-1; deriv>0; deriv--) {
      this->s[0].x[deriv-1] += _dt * (251.0*this->s[0].x[deriv] + 646.0*this->s[1].x[deriv] - 264.0*this->s[2].x[deriv] + 106.0*this->s[3].x[deriv] - 19.0*this->s[4].x[deriv]) / 720.0;
    }

    // perform forward integration: AB5 for velocity, AM5 for position
    //this->s[0].x[1] = this->s[1].x[1] +  _dt * (1901.0*this->s[1].x[2] - 2774.0*this->s[2].x[2] + 2616.0*this->s[3].x[2] - 1274.0*this->s[4].x[2] + 251.0*this->s[5].x[2]) / 720.0;
    //this->s[0].x[0] = this->s[1].x[0] +  _dt * (251.0*this->s[0].x[1] + 646.0*this->s[1].x[1] - 264.0*this->s[2].x[1] + 106.0*this->s[3].x[1] - 19.0*this->s[4].x[1]) / 720.0;

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * Standard Verlet (non-velocity) integrator
 */
template <class T>
class Verlet : public MultistepIntegrator<T> {
public:
  Verlet (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(2, _system, _level, _dt)
  {
    // prevent construction on a Velocity system?
    //assert(_system.hasAccel() && "Verlet cannot accept a VelocitySystem");
  }
  
  void stepForward (const double _dt) {
    assert(this->g.hasAccel() && "Verlet cannot integrate a VelocitySystem");

    // ask the system to find its new highest-level derivative
    //this->s[0].x[2] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration
    this->s[0].x[0] = 2.0*this->s[1].x[0] - this->s[2].x[0] +  _dt*_dt*this->s[1].x[2];

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * New integrator - uses more history, but still no velocities
 *
 * This is effectively performing Richardson extrapolation on the standard Verlet!
 */
template <class T>
class RichardsonVerlet : public MultistepIntegrator<T> {
public:
  RichardsonVerlet (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(4, _system, _level, _dt)
  {
    // prevent construction on a Velocity system?
    //assert(_system.hasAccel() && "RichardsonVerlet cannot accept a VelocitySystem");
  }

  void stepForward (const double _dt) {
    assert(this->g.hasAccel() && "RichardsonVerlet cannot integrate a VelocitySystem");

    // ask the system to find its new highest-level derivative
    //this->s[0].x[2] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration
    // note precision problems that will arise when positions are subtracted!
    // maybe the velocity version will be more precise?
    // just adding the parentheses increases accuracy!
    this->s[0].x[0] = this->s[3].x[0] + ((this->s[1].x[0] - this->s[4].x[0])
             + 0.25*_dt*_dt* ( 5.0*this->s[1].x[2]
                              +2.0*this->s[2].x[2]
                              +5.0*this->s[3].x[2]));

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * Hamming - pg 416
 *
 * A predictor using only accelerations and positions
 */
template <class T>
class Hamming416 : public MultistepIntegrator<T> {
public:
  Hamming416 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(4, _system, _level, _dt)
  {
    // prevent construction on a Velocity system?
    //assert(_system.hasAccel() && "Hamming416 cannot accept a VelocitySystem");
  }
  
  void stepForward (const double _dt) {
    assert(this->g.hasAccel() && "Hamming416 cannot integrate a VelocitySystem");

    // ask the system to find its new highest-level derivative
    //this->s[0].x[2] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead(this->g.getNumDerivs(), this->s[0].level, this->s[0].step++);
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration
    this->s[0].x[0] = 2.0*this->s[2].x[0] - this->s[4].x[0]
             + (4.0*_dt*_dt/3.0) * ( this->s[1].x[2]
                                    +this->s[2].x[2]
                                    +this->s[3].x[2]);

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * Hamming - pg 418
 *
 * A predictor using accelerations, velocities, and positions
 */
template <class T>
class Hamming418 : public MultistepIntegrator<T> {
public:
  Hamming418 (DynamicalSystem<T>& _system, const int32_t _level, const double _dt) :
    MultistepIntegrator<T>(2, _system, _level, _dt)
  {
    // prevent construction on a Velocity system?
    //assert(_system.hasAccel() && "Hamming418 cannot accept a VelocitySystem");
  }
  
  void stepForward (const double _dt) {
    assert(this->g.hasAccel() && "Hamming418 cannot integrate a VelocitySystem");

    // ask the system to find its new highest-level derivative
    //this->s[0].x[2] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0], this->getTime());

    // add a new one to the head
    DynamicState<T> newHead(this->g.getNumDerivs(), this->s[0].level, this->s[0].step++);
    newHead.time += _dt;
    this->s.insert(this->s.begin(), newHead);

    // perform forward integration
    // first, use Adams Moulton to find the new velocity
    // NEED SOMETHING BETTER HERE!
    //this->s[1].x[1] = this->s[2].x[1]
    //         + (_dt/2.0)      * (     this->s[1].x[2] +     this->s[2].x[2]);
    // then use the Hamming method to find the position
    this->s[0].x[0] = this->s[1].x[0]
             + (_dt/2.0)      * (    -this->s[1].x[1] + 3.0*this->s[2].x[1])
             + (_dt*_dt/12.0) * (17.0*this->s[1].x[2] + 7.0*this->s[2].x[2]);
    // this is horrible!
    this->s[0].x[1] = (384.0/_dt)    * (     this->s[2].x[0] - this->s[1].x[0])
             + (1.0)          * (312.*this->s[1].x[1] + 73.*this->s[2].x[1])
             + (-1.0*_dt)     * (110.*this->s[1].x[2] + 8.0*this->s[2].x[2]);

    // get rid of oldest
    this->s.pop_back();
  }
};


/*
 * Try Bulirsch-Stoer algorithm
 *
 * https://en.wikipedia.org/wiki/Bulirsch-Stoer_algorithm
*/

