#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
 * Verlet - a program to test forward advection schemes
 *
 * Copyright 2016 Mark J. Stock, markjstock@gmail.com
 *
 * v9: Dynamical state now has vector of derivatives
 *
 * Compile and run with:
 * g++ -o runmultistep -Ofast -std=c++11 -I/usr/include/eigen3 multistep.cpp && ./runmultistep
 */


/*
 * A class to contain one state (one set of changing variables, and all of their derivatives)
 */
class DynamicState {
public:
  // delegating constructor chain
  DynamicState() :
    DynamicState(1)
  {}

  DynamicState(const int _highestDeriv) :
    DynamicState(_highestDeriv, 0, 0)
  {}

  DynamicState(const int _highestDeriv, const int _level, const int _step) :
    x(1+_highestDeriv, ArrayXd()),
    level(_level), step(_step)
  {
    // do not initialize the arrays now, wait for later
  };

  // copy constructor works just like it's supposed to

  // like a copy constructor, but advances step and zeros highest derivatives
  DynamicState stepHelper() {
    // make the new object here
    DynamicState next(x.size(), level, step);
    // copy all derivatives except for the highest
    for (int i=0; i<x.size(); i++) {
      next.x[i] = x[i];
    }
    // highest derivative array is left as zeros
    next.x[x.size()] = ArrayXd::Zero(x[0].size());
    // finally, advance the step by one
    next.step++;
    return next;
  }

  double getFloatStep () {
    return step * pow(2.0, level);
  }

  ArrayXd getPos(void) {
    //cout << "DynamicState::getPos " << x.size() << endl;
    //cout << "DynamicState::getPos " << x.size() << "  " << x[0].size() << endl;
    //cout << "DynamicState::getPos " << x.size() << "  " << x[0].size() << "  " << x[0].segment(0,4).transpose() << endl;
    try {
      return x[0];
    } catch (exception& e) {
      cout << "Standard exception: " << e.what() << endl;
    }
    //if (x.size() > 0) return x[0];
    //else return ArrayXd::Zero(1);
  }

  ArrayXd getVel(void) {
    //if (x.size() > 1) return x[1];
    //else return ArrayXd::Zero(1);
    try {
      return x[1];
    } catch (exception& e) {
      cout << "Standard exception: " << e.what() << endl;
    }
  }

  ArrayXd getAcc(void) {
    //if (x.size() > 2) return x[2];
    //else return ArrayXd::Zero(1);
    try {
      return x[2];
    } catch (exception& e) {
      cout << "Standard exception: " << e.what() << endl;
    }
  }

  // level 0 is at base dt, level 1 is at 2*dt, level -1 at 0.5*dt
  int level;
  // step is time step within that level
  int step;
  // store a value and an arbitrary number of derivatives
  // x[0] is position, x[1] velocity, x[2] acceleration, x[3] jerk, etc.
  vector<ArrayXd> x;
};


/*
 * Must make a virtual superclass for integrable systems
 */
//class SineWave : public AccelerationSystem {
//class NBodyGrav : public AccelerationSystem {
//class NBodyVort : public VelocitySystem {

/*
 * DynamicalSystem - generalized class for dynamical systems
 */
class DynamicalSystem {
public:
  DynamicalSystem(const int _num) :
    numVars(_num)
  {}

  // number of derivatives in the dynamic state
  virtual int getNumDerivs(void) = 0;

  // return the initial state
  virtual DynamicState getInit(void) = 0;

  virtual bool hasAccel(void) = 0;
  virtual ArrayXd getHighestDeriv(const ArrayXd pos) = 0;

protected:
  // number of equations in the system
  int numVars;
};


/*
 * VelocitySystem - a dynamic system driven by velocities
 *                  like vortex methods; states have x, x'
 */
class VelocitySystem : public DynamicalSystem {
public:
  VelocitySystem(const int _num) :
    DynamicalSystem(_num),
    ic(DynamicState(1))
  {
    ic.step = 0;
  }

  // number of derivatives in the dynamic state
  int getNumDerivs(void) { return 1; }

  // return the initial state
  DynamicState getInit() { return ic; }

  bool hasAccel(void) { return false; }
  //virtual ArrayXd getHighestDeriv(const ArrayXd pos) = 0;

protected:
  // must define and store initial state
  DynamicState ic;
};


/*
 * AccelerationSystem - a dynamic system driven by accelerations
 *                      like gravitation; states have x, x', x"
 */
class AccelerationSystem : public DynamicalSystem {
public:
  AccelerationSystem(const int _num) :
    DynamicalSystem(_num),
    ic(DynamicState(2))
  {
    ic.step = 0;
  }

  // number of derivatives in the dynamic state
  int getNumDerivs(void) { return 2; }

  // return the initial state
  DynamicState getInit() { return ic; }

  bool hasAccel(void) { return true; }
  //virtual ArrayXd getHighestDeriv(const ArrayXd pos) = 0;

protected:
  // must define and store initial state
  DynamicState ic;
};


/*
 * A 2D vortex system with constant strengths, radii
 */
class NBodyVort2D : public VelocitySystem {
public:
  NBodyVort2D(const int _num) :
    num(_num),
    VelocitySystem(2*_num)
  {
    // num is how many bodies

    // we store the unchanging properties here
    circ = ArrayXd::Random(num) / (2.0*num);
    radiusSquared = 0.2 + 0.1 * ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();

    // store initial conditions in case we want to reuse this
    ic.x[0] = 10.0 * ArrayXd::Random(numVars);
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  ArrayXd getHighestDeriv(const ArrayXd pos) {
    ArrayXd newVels = ArrayXd::Zero(numVars);

    // evaluate using ispc-compiled subroutine
    // pos, circ, radiusSquared --> newVels

    // evaluate locally
    if (true) {
      for (int i=0; i<num; ++i) {
        // new velocities on particle i
        Vector2d thisVel(0.0, 0.0);
        for (int j=0; j<num; ++j) {
          // 20 flops
          // the influence of particle j
          Vector2d dx = pos.segment(2*j,2) - pos.segment(2*i,2);
          double invdist = 1.0/(dx.norm()+radiusSquared(j));
          double factor = circ(j) * invdist * invdist;
          thisVel[0] -= dx[1] * factor;
          thisVel[1] += dx[0] * factor;
        }
        newVels.segment(2*i,2) = thisVel;
      }
    }
    return newVels;
  }

private:
  // number of bodies
  int num;
  // these values are constant in time and unique to this system
  ArrayXd circ;
  ArrayXd radiusSquared;
};

/*
 * A gravitational n-body system
 */
class NBodyGrav : public AccelerationSystem {
public:
  NBodyGrav(const int _num) :
    num(_num),
    AccelerationSystem(3*_num)
  {
    // num is how many bodies

    // we store the unchanging properties here
    mass = (2.0 + ArrayXd::Random(num)) / (2.0*num);
    radiusSquared = 0.2 + 0.1 * ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();

    // store initial conditions in case we want to reuse this
    ic.x[0] = 10.0 * ArrayXd::Random(numVars);
    ic.x[1] = 1.0 * ArrayXd::Random(numVars);
    //cout << "NBodyGrav::NBodyGrav " << ic.x[0].segment(0,4).transpose() << endl;
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  ArrayXd getHighestDeriv(const ArrayXd pos) {
    ArrayXd newVal = ArrayXd::Zero(numVars);

    // evaluate using ispc-compiled subroutine
    // pos, mass, radiusSquared --> newVal

    // evaluate locally
    if (true) {
      for (int i=0; i<num; ++i) {
        // new accelerations on particle i
        Vector3d newAcc(0.0, 0.0, 0.0);
        for (int j=0; j<num; ++j) {
          // 20 flops
          // the influence of particle j
          Vector3d dx = pos.segment(3*j,3) - pos.segment(3*i,3);
          double invdist = 1.0/(dx.norm()+radiusSquared(j));
          newAcc += dx * (mass(j) * invdist * invdist * invdist);
        }
        newVal.segment(3*i,3) = newAcc;
      }
    }
    return newVal;
  }

private:
  // number of bodies
  int num;
  // these values are constant in time and unique to this system
  ArrayXd mass;
  ArrayXd radiusSquared;
};


/*
 * General forward integrator class
 */
class ForwardIntegrator {
public:
  ForwardIntegrator (DynamicalSystem& _system, const int _nstates, const int _level) :
    g(_system),
    s(_nstates, DynamicState(_system.getNumDerivs(), _level, 0))
  {
    // set the zero state
    s[0] = g.getInit();
  }

  // derived classes must define this method
  virtual void stepForward (const double _dt) = 0;

  // building block of many algorithms --- will only modify start if we solve for highest derivs
  DynamicState EulerStep(DynamicState initial, DynamicState derivs, const double _dt) {

    // create the state with copies of the lower derivatives
    DynamicState next = initial.stepHelper();

    // find the highest derivative - nope
    //initial.x[g.getNumDerivs()] = g.getHighestDeriv(start.x[0]);

    // step all lower derivatives forward
    for (int nd=0; nd<g.getNumDerivs(); nd++) {
      double factor = 1.0;
      for (int d=nd; d>=0; d--) {
        factor *= _dt / (nd-d+1);
        next.x[d] += factor * derivs.x[nd+1];
      }
    }

    // find the highest derivative - nope, again
    //next.x[g.getNumDerivs()] = g.getHighestDeriv(next.x[0]);

    // return the state
    return next;
  }

  DynamicState EulerStep(DynamicState start, const double _dt) {
    // call the general routine
    return EulerStep(start, start, _dt);
  }

  ArrayXd getPosition () {
    return s[0].getPos();
  }

  ArrayXd getVelocity () {
    return s[0].getVel();
  }

  ArrayXd getDeriv (const int deriv) {
    try {
      return s[0].x[deriv];
    } catch (exception& e) {
      cout << "Standard exception: " << e.what() << endl;
    }
  }

  double getError (const ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    return( temp.matrix().norm() / sqrt(temp.size()) );
  }

protected:
  DynamicalSystem& g;
  vector<DynamicState> s;
};


/*
 * Multi-stage forward integrator class (Euler, Runge-Kutta)
 */
class MultistageIntegrator : public ForwardIntegrator {
public:
  MultistageIntegrator (DynamicalSystem& _system, const int _level) :
    ForwardIntegrator(_system, 1, _level)
  {
    // zero state set in parent constructor
  }

protected:
  // save stages here?
  //vector<DynamicState> stage;
};


/*
 * Control class for Richardson Extrapolation and globally adaptive time stepping
 */
class RichardsonEuler {
public:
  RichardsonEuler (DynamicalSystem& _system) {

  }

private:
  //Euler
};


/*
 * Your basic 1-step Euler integrator
 */
class Euler : public MultistageIntegrator {
public:
  Euler (DynamicalSystem& _system, const int _level) :
    MultistageIntegrator(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    //int numDeriv = g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    s[0].x[g.getNumDerivs()] = g.getHighestDeriv(s[0].x[0]);

    // from that state, project forward
    DynamicState newHead = EulerStep(s[0], _dt);

    // add the new state to the head
    s.insert(s.begin(), newHead);

    // perform forward integration
    // s[1].x[0]  means  state t-1, 0th derivative

    /*
    for (int deriv=0; deriv<numDeriv; deriv++) {
      //s[0].x[deriv] = s[1].x[deriv];
      double factor = 1.0;
      for (int d=deriv; d>=0; d--) {
        factor *= _dt / (deriv-d+1);
        s[0].x[d] += factor * s[1].x[deriv+1];
      }
    }
    */

    /*
    s[0].x[0] = s[1].x[0];
    s[0].x[0] += _dt*s[1].x[1];
    if (g.hasAccel()) {
      s[0].x[1] = s[1].x[1];
      s[0].x[1] += _dt    *s[1].x[2];
      s[0].x[0] += _dt*_dt*s[1].x[2] / 2.0;
    }
    if (g.hasJerk()) {
      s[0].x[2] = s[1].x[2];
      s[0].x[2] += _dt        *s[1].x[3];
      s[0].x[1] += _dt*_dt    *s[1].x[3] / 2.0;
      s[0].x[0] += _dt*_dt*_dt*s[1].x[3] / 6.0;
    }
    */

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Runge-Kutta 2nd order
 */
class RK2 : public MultistageIntegrator {
public:
  RK2 (DynamicalSystem& _system, const int _level) :
    MultistageIntegrator(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // Generalizable to Heun's or Ralson's Methods
  void stepForward (const double _dt) {

    // alpha:error at 800 steps are  0.5:1.77703e-05  1.0:1.78663e-05  2/3:1.77534e-05
    const double alpha = 2.0/3.0;

    // ask the system to find its new highest-level derivative
    int numDeriv = g.getNumDerivs();
    s[0].x[numDeriv] = g.getHighestDeriv(s[0].x[0]);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    //const double hdt = alpha*_dt;
    //DynamicState stage2(numDeriv,0,0);

    // from that state, project forward
    DynamicState stage2 = EulerStep(s[0], alpha*_dt);
    stage2.x[numDeriv] = g.getHighestDeriv(stage2.x[0]);

    /*
    // if vel:
    stage2.x[0] = s[0].x[0] + hdt*s[0].x[1];
    stage2.x[1] = g.getHighestDeriv(stage2.x[0]);
    */

    // if accel:
    //stage2.x[0] = s[0].x[0] + hdt*s[0].x[1];
    // adding the last term halves the error!
    //stage2.x[0] = s[0].x[0] + hdt*s[0].x[1] + hdt*hdt*s[0].x[2]/2.0;
    //stage2.x[1] = s[0].x[1] + hdt*s[0].x[2];
    //stage2.x[2] = g.getHighestDeriv(stage2.x[0]);

    /*
    // if jerk:
    stage2.x[0] = s[0].x[0] + hdt*s[0].x[1] + hdt*hdt*s[0].x[2]/2.0 + hdt*hdt*hdt*s[0].x[3]/6.0;
    stage2.x[1] = s[0].x[1] + hdt*s[0].x[2] + hdt*hdt*s[0].x[3]/2.0;
    stage2.x[2] = s[0].x[2] + hdt*s[0].x[3];
    stage2.x[3] = g.getHighestDeriv(stage2.x[0]);
    */

    // add a new state to the head
    //DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    DynamicState newHead = s[0].stepHelper();

    // position updates via weighted average velocity
    const double oo2a = 1.0 / (2.0*alpha);
    for (int d=0; d<numDeriv; d++) {
      newHead.x[d] += _dt * ((1.0-oo2a)*s[0].x[d+1] + oo2a*stage2.x[d+1]);
      //s[0].x[d] = s[1].x[d] + _dt * ((1.0-oo2a)*s[1].x[d+1] + oo2a*stage2.x[d+1]);
    }
    //s[0].x[0] = s[1].x[0] + _dt * ((1.0-oo2a)*s[1].x[1] + oo2a*stage2.x[1]);
    //s[0].x[1] = s[1].x[1] + _dt * ((1.0-oo2a)*s[1].x[2] + oo2a*stage2.x[2]);

    // add a new state to the head
    s.insert(s.begin(), newHead);

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order
 */
class RK4 : public MultistageIntegrator {
public:
  RK4 (DynamicalSystem& _system, const int _level) :
    MultistageIntegrator(_system, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // Most write-ups of this are incorrect! Does nobody edit their shit?
  // Note that this could be improved using the 3/8 rule, see wikipedia
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    //cout << "in RK4::stepForward " << s[0].x[0].segment(0,4).transpose() << endl;

    // solve for top derivative at current state
    const int nd = g.getNumDerivs();
    s[0].x[nd] = g.getHighestDeriv(s[0].x[0]);

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
    DynamicState stage2 = s[0].stepHelper();
    for (int d=0; d<nd; d++) stage2.x[d] += hdt*s[0].x[d+1];
    //for (int d=0; d<nd; d++) stage2.x[d] = s[0].x[d] + hdt*s[0].x[d+1];
    //stage2.x[0] = s[0].x[0] + hdt*s[0].x[1];
    //stage2.x[1] = s[0].x[1] + hdt*s[0].x[2];
    stage2.x[nd] = g.getHighestDeriv(stage2.x[0]);

    // third step: project forward a half step from initial using the new acceleration
    DynamicState stage3 = s[0].stepHelper();
    for (int d=0; d<nd; d++) stage3.x[d] += hdt*stage2.x[d+1];
    //DynamicState stage3(nd,0,0);
    //stage3.x[0] = s[0].x[0] + hdt*stage2.x[1];
    //stage3.x[1] = s[0].x[1] + hdt*stage2.x[2];
    stage3.x[nd] = g.getHighestDeriv(stage3.x[0]);

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState stage4 = s[0].stepHelper();
    for (int d=0; d<nd; d++) stage4.x[d] += _dt*stage3.x[d+1];
    //DynamicState stage4(nd,0,0);
    //stage4.x[0] = s[0].x[0] + _dt*stage3.x[1];
    //stage4.x[1] = s[0].x[1] + _dt*stage3.x[2];
    stage4.x[nd] = g.getHighestDeriv(stage4.x[0]);

    // add a new state to the head
    //DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    DynamicState newHead = s[0].stepHelper();

    for (int d=0; d<nd; d++) {
      newHead.x[d] += _dt * (s[0].x[d+1] + 2.0*stage2.x[d+1] + 2.0*stage3.x[d+1] + stage4.x[d+1]) / 6.0;
    }

    // position updates via weighted average velocity
    //s[0].x[0] = s[1].x[0] + _dt * (s[1].x[1] + 2.0*stage2.x[1] + 2.0*stage3.x[1] + stage4.x[1]) / 6.0;
    // velocity updates via weighted average acceleration
    //s[0].x[1] = s[1].x[1] + _dt * (s[1].x[2] + 2.0*stage2.x[2] + 2.0*stage3.x[2] + stage4.x[2]) / 6.0;

    // add a new state to the head
    s.insert(s.begin(), newHead);

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Multi-step forward integrator class (saves previous solutions)
 */
class MultistepIntegrator : public ForwardIntegrator {
public:
  MultistepIntegrator (const int _nsteps, DynamicalSystem& _system, const int _level, const double _dt) :
    ForwardIntegrator(_system, _nsteps, _level)
  {
    // zero state set in parent constructor
    // set the previous states here
    RK4 r(g,0);
    //cout << "in MultistepIntegrator::MultistepIntegrator " << r.getPosition().segment(0,4).transpose() << endl;
    for (int istep=1; istep<_nsteps; ++istep) {
      for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
      s[istep].step = -istep;
      for (int deriv=0; deriv<g.getNumDerivs(); deriv++) {
        s[istep].x[deriv] = r.getDeriv(deriv);
      }
      s[istep].x[g.getNumDerivs()] = g.getHighestDeriv(s[istep].x[0]);
    }
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton
 */
class AB2 : public MultistepIntegrator {
public:
  AB2 (DynamicalSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(2, _system, _level, _dt)
  {}
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    // assert that we have at least two states in the history
    //static_assert(s.size() >= 2, "State vector does not have enough entries");
    // assert that this dt matches the previous dt
    const int numDerivs = g.getNumDerivs();

    // ask the system to find its new highest-level derivative
    s[0].x[numDerivs] = g.getHighestDeriv(s[0].x[0]);

    // add a new state to the head
    DynamicState newHead = s[0].stepHelper();
    s.insert(s.begin(), newHead);

    // perform forward integration: AB2 for first one down, AM2 for all others
    s[0].x[numDerivs-1] = s[1].x[numDerivs-1] + 0.5 * _dt * (3.0*s[1].x[numDerivs] - s[2].x[numDerivs]);
    for (int deriv=numDerivs-1; deriv>0; deriv--) {
      s[0].x[deriv-1] += 0.5 * _dt * (s[0].x[deriv] + s[1].x[deriv]);
    }

    // perform forward integration: AB2 for velocity, AM2 for position
    //s[0].x[1] = s[1].x[1] + 0.5 * _dt * (3.0*s[1].x[2] - s[2].x[2]);
    //s[0].x[0] = s[1].x[0] + 0.5 * _dt * (s[0].x[1] + s[1].x[1]);

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 4th order
 */
class AB4 : public MultistepIntegrator {
public:
  AB4 (DynamicalSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(4, _system, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    const int nd = g.getNumDerivs();

    s[0].x[nd] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead = s[0].stepHelper();
    s.insert(s.begin(), newHead);

    // perform forward integration: AB4 for first, AM4 for all lower-order derivatives
    s[0].x[nd-1] += _dt * (55.0*s[1].x[nd] - 59.0*s[2].x[nd] + 37.0*s[3].x[nd] - 9.0*s[4].x[nd]) / 24.0;
    for (int deriv=nd-1; deriv>0; deriv--) {
      s[0].x[deriv-1] += _dt * (9.0*s[0].x[deriv] + 19.0*s[1].x[deriv] - 5.0*s[2].x[deriv] + s[3].x[deriv]) / 24.0;
    }

    // perform forward integration: AB4 for velocity, AM4 for position
    //s[0].x[1] = s[1].x[1] +  _dt * (55.0*s[1].x[2] - 59.0*s[2].x[2] + 37.0*s[3].x[2] - 9.0*s[4].x[2]) / 24.0;
    //s[0].x[0] = s[1].x[0] +  _dt * (9.0*s[0].x[1] + 19.0*s[1].x[1] - 5.0*s[2].x[1] + s[3].x[1]) / 24.0;

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 5th order
 */
class AB5 : public MultistepIntegrator {
public:
  AB5 (DynamicalSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(5, _system, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    const int nd = g.getNumDerivs();

    s[0].x[nd] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead = s[0].stepHelper();
    s.insert(s.begin(), newHead);

    // perform forward integration: AB5 for first, AM5 for all lower-order derivatives
    s[0].x[nd-1] += _dt * (1901.0*s[1].x[nd] - 2774.0*s[2].x[nd] + 2616.0*s[3].x[nd] - 1274.0*s[4].x[nd] + 251.0*s[5].x[nd]) / 720.0;
    for (int deriv=nd-1; deriv>0; deriv--) {
      s[0].x[deriv-1] += _dt * (251.0*s[0].x[deriv] + 646.0*s[1].x[deriv] - 264.0*s[2].x[deriv] + 106.0*s[3].x[deriv] - 19.0*s[4].x[deriv]) / 720.0;
    }

    // perform forward integration: AB5 for velocity, AM5 for position
    //s[0].x[1] = s[1].x[1] +  _dt * (1901.0*s[1].x[2] - 2774.0*s[2].x[2] + 2616.0*s[3].x[2] - 1274.0*s[4].x[2] + 251.0*s[5].x[2]) / 720.0;
    //s[0].x[0] = s[1].x[0] +  _dt * (251.0*s[0].x[1] + 646.0*s[1].x[1] - 264.0*s[2].x[1] + 106.0*s[3].x[1] - 19.0*s[4].x[1]) / 720.0;

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Standard Verlet (non-velocity) integrator
 */
class Verlet : public MultistepIntegrator {
public:
  Verlet (AccelerationSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(2, _system, _level, _dt)
  {}
  
  void stepForward (const double _dt) {
    s[0].x[2] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    s[0].x[0] = 2.0*s[1].x[0] - s[2].x[0] +  _dt*_dt*s[1].x[2];

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * New integrator - uses more history, but still no velocities
 *
 * This is effectively performing Richardson extrapolation on the standard Verlet!
 */
class RichardsonVerlet : public MultistepIntegrator {
public:
  RichardsonVerlet (AccelerationSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(4, _system, _level, _dt)
  {}

  void stepForward (const double _dt) {
    s[0].x[2] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    // note precision problems that will arise when positions are subtracted!
    // maybe the velocity version will be more precise?
    // just adding the parentheses increases accuracy!
    s[0].x[0] = s[3].x[0] + ((s[1].x[0] - s[4].x[0])
             + 0.25*_dt*_dt* ( 5.0*s[1].x[2]
                              +2.0*s[2].x[2]
                              +5.0*s[3].x[2]));

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Hamming - pg 416
 *
 * A predictor using only accelerations and positions
 */
class Hamming416 : public MultistepIntegrator {
public:
  Hamming416 (AccelerationSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(4, _system, _level, _dt)
  {}
  
  void stepForward (const double _dt) {
    s[0].x[2] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    s[0].x[0] = 2.0*s[2].x[0] - s[4].x[0]
             + (4.0*_dt*_dt/3.0) * ( s[1].x[2]
                                    +s[2].x[2]
                                    +s[3].x[2]);

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Hamming - pg 418
 *
 * A predictor using accelerations, velocities, and positions
 */
class Hamming418 : public MultistepIntegrator {
public:
  Hamming418 (AccelerationSystem& _system, const int _level, const double _dt) :
    MultistepIntegrator(2, _system, _level, _dt)
  {}
  
  void stepForward (const double _dt) {
    s[0].x[2] = g.getHighestDeriv(s[0].x[0]);

    // add a new one to the head
    DynamicState newHead(g.getNumDerivs(), s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    // first, use Adams Moulton to find the new velocity
    // NEED SOMETHING BETTER HERE!
    //s[1].x[1] = s[2].x[1]
    //         + (_dt/2.0)      * (     s[1].x[2] +     s[2].x[2]);
    // then use the Hamming method to find the position
    s[0].x[0] = s[1].x[0]
             + (_dt/2.0)      * (    -s[1].x[1] + 3.0*s[2].x[1])
             + (_dt*_dt/12.0) * (17.0*s[1].x[2] + 7.0*s[2].x[2]);
    // this is horrible!
    s[0].x[1] = (384.0/_dt)    * (     s[2].x[0] - s[1].x[0])
             + (1.0)          * (312.*s[1].x[1] + 73.*s[2].x[1])
             + (-1.0*_dt)     * (110.*s[1].x[2] + 8.0*s[2].x[2]);

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Try Bulirsch-Stoer algorithm
 *
 * https://en.wikipedia.org/wiki/Bulirsch-Stoer_algorithm
*/


// Perform Richardson Extrapolation by creating multiple temporal resolution levels automatically?


// Create a system and an integrator
int main () {

  // iterate a gravitational n-body system for a few steps
  NBodyGrav s(100);
  //NBodyVort2D vort(100);
  //Euler ev(vort,0);

  // find the "exact" solution for AB4 - use this as the exact solution for everyone
  double time = 0.0;
  int maxSteps = 10000;
  double dt = 10.0 / maxSteps;
  AB5 exact(s,0,dt);
  cout << "'Exact' solution is from running " << maxSteps << " steps of AB5 at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    //ev.stepForward(dt);
    time += dt;
  }

  cout << "steps\t";
  cout << "Euler\t\t";
  cout << "RK2\t\t";
  cout << "AB2\t\t";
  cout << "Verlet\t\t";
  cout << "RK4\t\t";
  cout << "AB4\t\t";
  cout << "Verlet4\t\t";
  cout << "Ham416\t\t";
  //cout << "Ham418\t\t";
  cout << "AB5";
  cout << endl;

  // integrate using the various methods
  for (maxSteps = 12; maxSteps < 15000; maxSteps *= 2) {
    dt = 10.0 / maxSteps;
    //cout << "Running " << maxSteps << " steps at dt= " << dt << endl;
    time = 0.0;

    // initialize integrators
    Euler e(s,0);
    RK2 r2(s,0);
    AB2 a(s,0,dt);
    Verlet ve(s,0,dt);
    RK4 r(s,0);
    AB4 ab(s,0,dt);
    RichardsonVerlet rv(s,0,dt);
    Hamming416 h6(s,0,dt);
    Hamming418 h8(s,0,dt);
    AB5 abb(s,0,dt);

    for (int i=0; i<maxSteps; ++i) {

      // take the forward step
      e.stepForward(dt);
      if (i%2 == 0) r2.stepForward(2.0*dt);
      a.stepForward(dt);
      ve.stepForward(dt);
      if (i%4 == 0) r.stepForward(4.0*dt);
      ab.stepForward(dt);
      rv.stepForward(dt);
      h6.stepForward(dt);
      h8.stepForward(dt);
      abb.stepForward(dt);
      time += dt;

    }

    ArrayXd eSolution = e.getPosition();
    //cout << "  Euler solution: " << eSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in Euler is " << exact.getError(eSolution) << endl;

    ArrayXd r2Solution = r2.getPosition();
    //cout << "  RK2 solution: " << r2Solution.segment(0,4).transpose() << endl;
    //cout << "  Error in RK2 is " << exact.getError(r2Solution) << endl;

    ArrayXd aSolution = a.getPosition();
    //cout << "  ABM2 solution: " << aSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in ABM2 is " << exact.getError(aSolution) << endl;

    ArrayXd vSolution = ve.getPosition();
    //cout << "  Verlet solution: " << vSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in Verlet is " << exact.getError(vSolution) << endl;

    ArrayXd rSolution = r.getPosition();
    //cout << "  RK4 solution: " << rSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in RK4 is " << exact.getError(rSolution) << endl;

    ArrayXd abSolution = ab.getPosition();
    //cout << "  ABM4 solution: " << abSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in ABM4 is " << exact.getError(abSolution) << endl;

    ArrayXd sSolution = rv.getPosition();
    //cout << "  RichardsonVerlet solution: " << sSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in RichardsonVerlet is " << exact.getError(sSolution) << endl;

    ArrayXd h6Solution = h6.getPosition();
    ArrayXd h8Solution = h8.getPosition();
    ArrayXd abbSolution = abb.getPosition();

    cout << maxSteps << "\t" << exact.getError(eSolution);
    cout << "\t" << exact.getError(r2Solution);
    cout << "\t" << exact.getError(aSolution);
    cout << "\t" << exact.getError(vSolution);
    cout << "\t" << exact.getError(rSolution);
    cout << "\t" << exact.getError(abSolution);
    cout << "\t" << exact.getError(sSolution);
    cout << "\t" << exact.getError(h6Solution);
    //cout << "\t" << exact.getError(h8Solution);
    cout << "\t" << exact.getError(abbSolution);
    cout << endl;
  }

  exit(0);
}
