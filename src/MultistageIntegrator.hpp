/*
 * MultistageIntegrator.hpp - Multi-stage forward integrator classes (Euler, Runge-Kutta)
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#pragma once

#include "DynamicState.hpp"
#include "DynamicalSystem.hpp"
#include "ForwardIntegrator.hpp"

#include <cstdint>
#include <iostream>
#include <vector>


/*
 * Multi-stage forward integrator class (Euler, Runge-Kutta)
 */
template <class T>
class MultistageIntegrator : public ForwardIntegrator<T> {
public:
  MultistageIntegrator (DynamicalSystem<T>& _system, const int32_t _level) :
    ForwardIntegrator<T>(_system, 1, _level)
  {
    // zero state set in parent constructor
  }

protected:
  // save stages in derived classes
};


/*
 * Your basic 1-step explicit Euler integrator
 *
 *  0  | 
 * ---------
 *     |  1
 */
template <class T>
class Euler : public MultistageIntegrator<T> {
public:
  Euler (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> newHead = this->CombinedStepEnhanced(this->s[0], _dt, 1, this->s[0], 1.0);

    // add the new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 2nd order - Heun's Method
 *
 *  0  | 
 *  1  |  1
 * -------------
 *     | 1/2 1/2
 */
template <class T>
class RK2Heun : public MultistageIntegrator<T> {
public:
  RK2Heun (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> stage2 = this->CombinedStepEnhanced(this->s[0], _dt, 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // or do it this way
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 2, this->s[0], 0.5, stage2, 0.5);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

/*
 * Runge-Kutta 2nd order - Ralston's
 *
 *  0  | 
 * 2/3 | 2/3
 * -------------
 *     | 1/4 3/4
 */
template <class T>
class RK2Ralston : public MultistageIntegrator<T> {
public:
  RK2Ralston (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> stage2 = this->CombinedStepEnhanced(this->s[0], _dt*2./3., 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // add it on this way
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 2, this->s[0], 0.25, stage2, 0.75);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

/*
 * Runge-Kutta 3rd order - Kutta's
 *
 *  0  | 
 * 1/2 | 1/2
 *  1  | -1   2 
 * ------------------
 *     | 1/6 2/3 1/6
 */
template <class T>
class RK3Kutta : public MultistageIntegrator<T> {
public:
  RK3Kutta (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> stage2 = this->CombinedStep(this->s[0], 0.5*_dt, 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, two derivs)
    // need to add -1*stage1 and 2*stage2 to get test positions: stage3
    DynamicState<T> stage3 = this->CombinedStep(this->s[0], _dt, 2, this->s[0], -1.0, stage2, 2.0);
    this->g.setHighestDeriv(stage3);

    // combine stages into new state
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 3, this->s[0], 1.0/6.0, stage2, 4.0/6.0, stage3, 1.0/6.0);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 3rd order - Ralston's method
 *
 *  0  | 
 * 1/2 | 1/2
 * 3/4 |  0  3/4
 * ------------------
 *     | 2/9 1/3 4/9
 */
template <class T>
class RK3Ralston : public MultistageIntegrator<T> {
public:
  RK3Ralston (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> stage2 = this->CombinedStep(this->s[0], 0.5*_dt, 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, new derivs)
    DynamicState<T> stage3 = this->CombinedStep(this->s[0], 0.75*_dt, 1, stage2, 1.0);
    this->g.setHighestDeriv(stage3);

    // combine stages into new state
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 3, this->s[0], 2./9., stage2, 1./3., stage3, 4./9.);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 3rd order - Heun's RK3
 *
 *  0  | 
 * 1/3 | 1/3
 * 2/3 |  0  2/3
 * ------------------
 *     | 1/4  0  3/4
 */
template <class T>
class RK3Heun : public MultistageIntegrator<T> {
public:
  RK3Heun (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> stage2 = this->CombinedStepEnhanced(this->s[0], _dt/3., 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, new derivs)
    DynamicState<T> stage3 = this->CombinedStep(this->s[0], _dt*2./3., 1, stage2, 1.0);
    this->g.setHighestDeriv(stage3);

    // combine stages into new state
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 2, this->s[0], 1./4., stage3, 3./4.);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order - classic
 *
 *  0  | 
 * 1/2 | 1/2
 * 1/2 |  0  1/2
 *  1  |  0   0   1
 * ----------------------
 *     | 1/6 1/3 1/3 1/6
 */
template <class T>
class RK4 : public MultistageIntegrator<T> {
public:
  RK4 (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // Most write-ups of this are incorrect! Does nobody edit their shit?
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    //std::cout << "in RK4::stepForward " << s[0].x[0].segment(0,4).transpose() << std::endl;

    // solve for top derivative at current state
    this->g.setHighestDeriv(this->s[0]);

    // second step: project forward a half step using that acceleration
    DynamicState<T> stage2 = this->CombinedStep(this->s[0], 0.5*_dt, 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // third step: project forward a half step from initial using the new acceleration
    DynamicState<T> stage3 = this->CombinedStep(this->s[0], 0.5*_dt, 1, stage2, 1.0);
    this->g.setHighestDeriv(stage3);

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState<T> stage4 = this->CombinedStep(this->s[0], _dt, 1, stage3, 1.0);
    this->g.setHighestDeriv(stage4);

    // new style
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 4, this->s[0], 1./6., stage2, 2./6., stage3, 2./6., stage4, 1./6.);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order - 3/8ths Rule - NOT DONE
 *
 *  0  | 
 * 1/3 | 1/3
 * 2/3 |-1/3  1
 *  1  |  1  -1   1
 * ----------------------
 *     | 1/8 3/8 3/8 1/8
 */
template <class T>
class RK4ter : public MultistageIntegrator<T> {
public:
  RK4ter (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
    //assert(false and "RK4ter incomplete!");
  }
  
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    //std::cout << "in RK4ter::stepForward " << s[0].x[0].segment(0,4).transpose() << std::endl;

    // solve for top derivative at current state
    this->g.setHighestDeriv(this->s[0]);

    // second step: project forward a half step using that acceleration
    DynamicState<T> stage2 = this->CombinedStep(this->s[0], _dt*1./3., 1, this->s[0], 1.0);
    this->g.setHighestDeriv(stage2);

    // third step: project forward a half step from initial using the new acceleration
    DynamicState<T> stage3 = this->CombinedStep(this->s[0], _dt*2./3., 2, this->s[0], -0.5, stage2, 1.5);
    this->g.setHighestDeriv(stage3);

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState<T> stage4 = this->CombinedStep(this->s[0], _dt, 3, this->s[0], 1.0, stage2, -1.0, stage3, 1.0);
    this->g.setHighestDeriv(stage4);

    // new style
    DynamicState<T> newHead = this->CombinedStep(this->s[0], _dt, 4, this->s[0], 1./8., stage2, 3./8., stage3, 3./8., stage4, 1./8.);
    newHead.time += _dt;

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

