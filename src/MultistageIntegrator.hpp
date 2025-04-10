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
 * Your basic 1-step Euler integrator
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
    //this->s[0].x[this->g.getNumDerivs()] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // from that state, project forward
    DynamicState<T> newHead = this->EulerStep(this->s[0], _dt);

    // add the new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 2nd order - Heun's Method
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
    //int32_t numDeriv = this->g.getNumDerivs();
    //this->s[0].x[numDeriv] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    //const double hdt = alpha*_dt;
    //DynamicState stage2(numDeriv,0,0);

    // from that state, project forward
    DynamicState<T> stage2 = this->EulerStep(this->s[0], _dt);
    //stage2.x[numDeriv] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+_dt);
    this->g.setHighestDeriv(stage2);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position updates via weighted average velocity
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (0.5*this->s[0].x[d+1] + 0.5*stage2.x[d+1]);
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

/*
 * Runge-Kutta 2nd order - Ralston's
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
    //int32_t numDeriv = this->g.getNumDerivs();
    //this->s[0].x[numDeriv] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // from that state, project forward
    DynamicState<T> stage2 = this->EulerStep(this->s[0], 2.0/3.0*_dt);
    //stage2.x[numDeriv] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+2.0/3.0*_dt);
    this->g.setHighestDeriv(stage2);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position updates via weighted average velocity
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (0.25*this->s[0].x[d+1] + 0.75*stage2.x[d+1]);
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

/*
 * Runge-Kutta 3rd order - Classic
 */
template <class T>
class RK3Kutta : public MultistageIntegrator<T> {
public:
  RK3Kutta (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
    assert(false and "RK3Kutta incomplete!");
  }
  
  void stepForward (const double _dt) {

    // ask the system to find its new highest-level derivative
    //const int32_t nd = this->g.getNumDerivs();
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // from that state, project forward
    DynamicState<T> stage2 = this->EulerStep(this->s[0], 0.5*_dt);
    //stage2.x[nd] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+0.5*_dt);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, new derivs)
    // NEED A SMARTER WAY TO DO THIS - I AM CONFUSED!
    // need to add -1*stage1 and 2*stage2 to get test positions: stage3
    DynamicState<T> stage3 = this->EulerStep(this->s[0], stage2, _dt);
    //stage3.x[nd] = this->g.getHighestDeriv(stage3.x[0], this->getTime()+_dt);
    this->g.setHighestDeriv(stage3);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position updates via weighted average velocity
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (this->s[0].x[d+1] + 4.0*stage2.x[d+1] + stage3.x[d+1]) / 6.0;
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 3rd order - Ralston's method
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

    // first step: set stage 1 to the last solution (now s[1])

    // from that state, project forward
    DynamicState<T> stage2 = this->EulerStep(this->s[0], 0.5*_dt);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, new derivs)
    DynamicState<T> stage3 = this->EulerStep(this->s[0], stage2, 0.75*_dt);
    this->g.setHighestDeriv(stage3);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position updates via weighted average velocity
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (2.0*this->s[0].x[d+1] + 3.0*stage2.x[d+1] + 4.0*stage3.x[d+1]) / 9.0;
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 3rd order - Heun's RK3
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
    //const int32_t nd = this->g.getNumDerivs();
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // from that state, project forward
    DynamicState<T> stage2 = this->EulerStep(this->s[0], _dt/3.0);
    //stage2.x[nd] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+_dt/3.0);
    this->g.setHighestDeriv(stage2);

    // and do it again (using initial positions, new derivs)
    DynamicState<T> stage3 = this->EulerStep(this->s[0], stage2, _dt*2.0/3.0);
    //stage3.x[nd] = this->g.getHighestDeriv(stage3.x[0], this->getTime()+_dt*2.0/3.0);
    this->g.setHighestDeriv(stage3);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position updates via weighted average velocity
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (this->s[0].x[d+1] + 3.0*stage3.x[d+1]) / 4.0;
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order - classic
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
    //const int32_t nd = this->g.getNumDerivs();
    //this->s[0].x[nd] = this->g.getHighestDeriv(this->s[0].x[0], this->getTime());
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])
    const double hdt = 0.5*_dt;

    // This new way of doing the calculation DOES NOT HELP!
    // second step: project forward a half step using that acceleration
    //DynamicState stage2 = EulerStep(s[0], hdt);
    //stage2.x[nd] = g.getHighestDeriv(stage2.x[0]);

    // third step: project forward a half step from initial using the new acceleration
    //DynamicState stage3 = EulerStep(s[0], stage2, hdt);
    //stage3.x[nd] = g.getHighestDeriv(stage3.x[0]);

    // fourth step: project forward a full step from initial using the newest acceleration
    //DynamicState stage4 = EulerStep(s[0], stage3, _dt);
    //stage4.x[nd] = g.getHighestDeriv(stage4.x[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    DynamicState<T> stage2 = this->s[0].stepHelper();
    stage2.time += hdt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage2.x[d] += hdt*this->s[0].x[d+1];
    //for (int32_t d=0; d<nd; d++) stage2.x[d] = s[0].x[d] + hdt*s[0].x[d+1];
    //stage2.x[0] = s[0].x[0] + hdt*s[0].x[1];
    //stage2.x[1] = s[0].x[1] + hdt*s[0].x[2];
    //stage2.x[nd] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+hdt);
    this->g.setHighestDeriv(stage2);

    // third step: project forward a half step from initial using the new acceleration
    DynamicState<T> stage3 = this->s[0].stepHelper();
    stage3.time += hdt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage3.x[d] += hdt*stage2.x[d+1];
    //DynamicState stage3(nd,0,0);
    //stage3.x[0] = s[0].x[0] + hdt*stage2.x[1];
    //stage3.x[1] = s[0].x[1] + hdt*stage2.x[2];
    //stage3.x[nd] = this->g.getHighestDeriv(stage3.x[0], this->getTime()+hdt);
    this->g.setHighestDeriv(stage3);

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState<T> stage4 = this->s[0].stepHelper();
    stage4.time += _dt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage4.x[d] += _dt*stage3.x[d+1];
    //DynamicState stage4(nd,0,0);
    //stage4.x[0] = s[0].x[0] + _dt*stage3.x[1];
    //stage4.x[1] = s[0].x[1] + _dt*stage3.x[2];
    //stage4.x[nd] = this->g.getHighestDeriv(stage4.x[0], this->getTime()+_dt);
    this->g.setHighestDeriv(stage4);

    // add a new state to the head
    //DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position, vel, etc. updates use weighted averages
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (this->s[0].x[d+1] + 2.0*stage2.x[d+1] + 2.0*stage3.x[d+1] + stage4.x[d+1]) / 6.0;
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order - 3/8ths Rule - NOT DONE
 */
template <class T>
class RK4ter : public MultistageIntegrator<T> {
public:
  RK4ter (DynamicalSystem<T>& _system, const int32_t _level) :
    MultistageIntegrator<T>(_system, _level)
  {
    // initial conditions set in parent constructor
    assert(false and "RK4ter incomplete!");
  }
  
  // Most write-ups of this are incorrect! Does nobody edit their shit?
  // Note that this could be improved using the 3/8 rule, see wikipedia
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    //std::cout << "in RK4ter::stepForward " << s[0].x[0].segment(0,4).transpose() << std::endl;

    // solve for top derivative at current state
    this->g.setHighestDeriv(this->s[0]);

    // first step: set stage 1 to the last solution (now s[1])
    const double hdt = 0.5*_dt;

    // This new way of doing the calculation DOES NOT HELP!
    // second step: project forward a half step using that acceleration
    //DynamicState stage2 = EulerStep(s[0], hdt);
    //stage2.x[nd] = g.getHighestDeriv(stage2.x[0]);

    // third step: project forward a half step from initial using the new acceleration
    //DynamicState stage3 = EulerStep(s[0], stage2, hdt);
    //stage3.x[nd] = g.getHighestDeriv(stage3.x[0]);

    // fourth step: project forward a full step from initial using the newest acceleration
    //DynamicState stage4 = EulerStep(s[0], stage3, _dt);
    //stage4.x[nd] = g.getHighestDeriv(stage4.x[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    DynamicState<T> stage2 = this->s[0].stepHelper();
    stage2.time += hdt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage2.x[d] += hdt*this->s[0].x[d+1];
    //for (int32_t d=0; d<nd; d++) stage2.x[d] = s[0].x[d] + hdt*s[0].x[d+1];
    //stage2.x[0] = s[0].x[0] + hdt*s[0].x[1];
    //stage2.x[1] = s[0].x[1] + hdt*s[0].x[2];
    //stage2.x[nd] = this->g.getHighestDeriv(stage2.x[0], this->getTime()+hdt);
    this->g.setHighestDeriv(stage2);

    // third step: project forward a half step from initial using the new acceleration
    DynamicState<T> stage3 = this->s[0].stepHelper();
    stage3.time += hdt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage3.x[d] += hdt*stage2.x[d+1];
    //DynamicState stage3(nd,0,0);
    //stage3.x[0] = s[0].x[0] + hdt*stage2.x[1];
    //stage3.x[1] = s[0].x[1] + hdt*stage2.x[2];
    //stage3.x[nd] = this->g.getHighestDeriv(stage3.x[0], this->getTime()+hdt);
    this->g.setHighestDeriv(stage3);

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState<T> stage4 = this->s[0].stepHelper();
    stage4.time += _dt;
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) stage4.x[d] += _dt*stage3.x[d+1];
    //DynamicState stage4(nd,0,0);
    //stage4.x[0] = s[0].x[0] + _dt*stage3.x[1];
    //stage4.x[1] = s[0].x[1] + _dt*stage3.x[2];
    //stage4.x[nd] = this->g.getHighestDeriv(stage4.x[0], this->getTime()+_dt);
    this->g.setHighestDeriv(stage4);

    // add a new state to the head
    DynamicState<T> newHead = this->s[0].stepHelper();
    newHead.time += _dt;

    // position, vel, etc. updates use weighted averages
    for (int32_t d=0; d<this->g.getNumDerivs(); d++) {
      newHead.x[d] += _dt * (this->s[0].x[d+1] + 2.0*stage2.x[d+1] + 2.0*stage3.x[d+1] + stage4.x[d+1]) / 6.0;
    }

    // add a new state to the head
    this->s.insert(this->s.begin(), newHead);

    // get rid of oldest state
    this->s.pop_back();
  }
};

