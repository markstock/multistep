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
 * v7: Refactoring integrators into classes
 *
 * Compile and run with:
 * g++ -o verlet -Ofast -std=c++11 -I/usr/include/eigen3 verlet7.cpp && ./verlet
 */


/*
 * A class to contain one state (one set of changing variables, and all of their derivatives)
 */
class DynamicState {
public:
  DynamicState(const int _level, const int _step) :
    level(_level), step(_step)
  {
    // do not initialize the arrays yet, wait for later
  };

  ~DynamicState() {};

  double getFloatStep () {
    return step * pow(2.0, level);
  }

//private:
  // level 0 is at base dt, level 1 is at 2*dt, level -1 at 0.5*dt
  int level;
  // step is time step within that level
  int step;
  // the three arrays: value, first and second derivatives
  ArrayXd pos;
  ArrayXd vel;
  ArrayXd acc;
  // ultimately, we want to store a value and an arbitrary number of derivatives
  // x[0] is position, x[1] velocity, x[2] acceleration, x[3] jerk, etc.
  //vector<DynamicState> x;
};


/*
 * Must make a virtual superclass for integrable systems
 */
//virtual class IntegrableSystem {
//class SineWave : public IntegrableSystem {
//class NBodyGrav : public IntegrableSystem {


/*
 * A gravitational n-body system
 */
class NBodyGrav {
public:
  NBodyGrav(const int _num) {
    // this is how many bodies
    num = _num;
    // we store the unchanging properties here
    mass = (2.0 + ArrayXd::Random(num)) / (2.0*num);
    radiusSquared = 0.2 + 0.1 * ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();
    // store initial conditions in case we want to reuse this
    initPos = 10.0 * ArrayXd::Random(3*num);
    initVel = 1.0 * ArrayXd::Random(3*num);
  };

  ~NBodyGrav() {};

  // initialize the system
  ArrayXd initPosit() {
    return initPos;
  }

  ArrayXd initVeloc() {
    return initVel;
  }

  // perform n-body acceleration calculation; uses position and mass and radius squared
  ArrayXd getAccel(const ArrayXd pos) {
    ArrayXd newVal = ArrayXd::Zero(3*num);

    // evaluate using ispc-compiled subroutine
    // pos, mass, radiusSquared --> newVal

    // evaluate locally
    if (true) {
      for (int i=0; i<num; ++i) {
        // new accelerations on particle i
        Vector3d temp = pos.segment(3*i,3);
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
  int num;
  // these values are constant in time
  ArrayXd mass;
  ArrayXd radiusSquared;
  ArrayXd initPos;
  ArrayXd initVel;
};


/*
 * General forward integrator class
 */
class ForwardIntegrator {
public:
  ForwardIntegrator (NBodyGrav& _gravSys, const int _nstates, const int _level) :
    g(_gravSys),
    s(_nstates, DynamicState(_level, 0))
  {
    // set the zero state
    s[0].step = 0;
    s[0].pos = g.initPosit();
    s[0].vel = g.initVeloc();
  }

  // derived classes must define this method
  virtual void stepForward (const double _dt) = 0;

  ArrayXd getPosition () {
    return s[0].pos;
  }

  ArrayXd getVelocity () {
    return s[0].vel;
  }

  double getError (const ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    return( temp.matrix().norm() / sqrt(temp.size()) );
  }

protected:
  NBodyGrav& g;
  vector<DynamicState> s;
};


/*
 * Multi-stage forward integrator class (Euler, Runge-Kutta)
 */
class MultistageIntegrator : public ForwardIntegrator {
public:
  MultistageIntegrator (NBodyGrav& _gravSys, const int _level) :
    ForwardIntegrator(_gravSys, 1, _level)
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
  RichardsonEuler (NBodyGrav& _gravSys) {

  }

private:
  //Euler
};


/*
 * Your basic 1-step Euler integrator
 */
class Euler : public MultistageIntegrator {
public:
  Euler (NBodyGrav& _gravSys, const int _level) :
    MultistageIntegrator(_gravSys, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    s[0].acc = g.getAccel(s[0].pos);

    // add a new state to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    s[0].pos = s[1].pos + _dt*s[1].vel + 0.5*_dt*_dt*s[1].acc;
    s[0].vel = s[1].vel + _dt*s[1].acc;

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Runge-Kutta 2nd order
 */
class RK2 : public MultistageIntegrator {
public:
  RK2 (NBodyGrav& _gravSys, const int _level) :
    MultistageIntegrator(_gravSys, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // Generalizable to Heun's or Ralson's Methods
  void stepForward (const double _dt) {
    // alpha:error at 800 steps are  0.5:1.77703e-05  1.0:1.78663e-05  2/3:1.77534e-05
    const double alpha = 2.0/3.0;

    // ask the system to find its new highest-level derivative
    s[0].acc = g.getAccel(s[0].pos);

    // add a new state to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    const double hdt = alpha*_dt;
    DynamicState stage2(0,0);
    stage2.pos = s[1].pos + hdt*s[1].vel;
    stage2.acc = g.getAccel(stage2.pos);
    stage2.vel = s[1].vel + hdt*s[1].acc;

    // position updates via weighted average velocity
    const double oo2a = 1.0 / (2.0*alpha);
    s[0].pos = s[1].pos + _dt * ((1.0-oo2a)*s[1].vel + oo2a*stage2.vel);
    // velocity updates via weighted average acceleration
    s[0].vel = s[1].vel + _dt * ((1.0-oo2a)*s[1].acc + oo2a*stage2.acc);

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Runge-Kutta 4th order
 */
class RK4 : public MultistageIntegrator {
public:
  RK4 (NBodyGrav& _gravSys, const int _level) :
    MultistageIntegrator(_gravSys, _level)
  {
    // initial conditions set in parent constructor
  }
  
  // Most write-ups of this are incorrect! Does nobody edit their shit?
  // Note that this could be improved using the 3/8 rule, see wikipedia
  void stepForward (const double _dt) {
    // ask the system to find its new highest-level derivative
    s[0].acc = g.getAccel(s[0].pos);

    // add a new state to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // first step: set stage 1 to the last solution (now s[1])

    // second step: project forward a half step using that acceleration
    const double hdt = 0.5*_dt;
    DynamicState stage2(0,0);
    stage2.pos = s[1].pos + hdt*s[1].vel;
    stage2.acc = g.getAccel(stage2.pos);
    stage2.vel = s[1].vel + hdt*s[1].acc;

    // third step: project forward a half step from initial using the new acceleration
    DynamicState stage3(0,0);
    stage3.pos = s[1].pos + hdt*stage2.vel;
    stage3.acc = g.getAccel(stage3.pos);
    stage3.vel = s[1].vel + hdt*stage2.acc;

    // fourth step: project forward a full step from initial using the newest acceleration
    DynamicState stage4(0,0);
    stage4.pos = s[1].pos + _dt*stage3.vel;
    stage4.acc = g.getAccel(stage4.pos);
    stage4.vel = s[1].vel + _dt*stage3.acc;

    // position updates via weighted average velocity
    s[0].pos = s[1].pos + _dt * (s[1].vel + 2.0*stage2.vel + 2.0*stage3.vel + stage4.vel) / 6.0;
    // velocity updates via weighted average acceleration
    s[0].vel = s[1].vel + _dt * (s[1].acc + 2.0*stage2.acc + 2.0*stage3.acc + stage4.acc) / 6.0;
  }
};


/*
 * Multi-step forward integrator class (saves previous solutions)
 */
class MultistepIntegrator : public ForwardIntegrator {
public:
  MultistepIntegrator (const int _nsteps, NBodyGrav& _gravSys, const int _level, const double _dt) :
    ForwardIntegrator(_gravSys, _nsteps, _level)
  {
    // zero state set in parent constructor
    // set the previous states here
    RK4 r(g,0);
    for (int istep=1; istep<_nsteps; ++istep) {
      for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
      s[istep].step = -istep;
      s[istep].pos = r.getPosition();
      s[istep].vel = r.getVelocity();
      s[istep].acc = g.getAccel(s[istep].pos);
    }
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton
 */
class AB2 : public MultistepIntegrator {
public:
  AB2 (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(2, _gravSys, _level, _dt)
  {}
  
  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    // assert that we have at least two states in the history
    //static_assert(s.size() >= 2, "State vector does not have enough entries");
    // assert that this dt matches the previous dt

    // ask the system to find its new highest-level derivative
    s[0].acc = g.getAccel(s[0].pos);

    // add a new state to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration: AB2 for velocity, AM2 for position
    s[0].vel = s[1].vel + 0.5 * _dt * (3.0*s[1].acc - s[2].acc);
    s[0].pos = s[1].pos + 0.5 * _dt * (s[0].vel + s[1].vel);

    // get rid of oldest state
    s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 4th order
 */
class AB4 : public MultistepIntegrator {
public:
  AB4 (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(4, _gravSys, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration: AB4 for velocity, AM4 for position
    s[0].vel = s[1].vel +  _dt * (55.0*s[1].acc - 59.0*s[2].acc + 37.0*s[3].acc - 9.0*s[4].acc) / 24.0;
    s[0].pos = s[1].pos +  _dt * (9.0*s[0].vel + 19.0*s[1].vel - 5.0*s[2].vel + s[3].vel) / 24.0;

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Adams-Bashforth plus Adams-Moulton, all 5th order
 */
class AB5 : public MultistepIntegrator {
public:
  AB5 (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(5, _gravSys, _level, _dt)
  {}

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration: AB5 for velocity, AM5 for position
    s[0].vel = s[1].vel +  _dt * (1901.0*s[1].acc - 2774.0*s[2].acc + 2616.0*s[3].acc - 1274.0*s[4].acc + 251.0*s[5].acc) / 720.0;
    s[0].pos = s[1].pos +  _dt * (251.0*s[0].vel + 646.0*s[1].vel - 264.0*s[2].vel + 106.0*s[3].vel - 19.0*s[4].vel) / 720.0;

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Standard Verlet (non-velocity) integrator
 */
class Verlet : public MultistepIntegrator {
public:
  Verlet (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(2, _gravSys, _level, _dt)
  {}
  
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    s[0].pos = 2.0*s[1].pos - s[2].pos +  _dt*_dt*s[1].acc;

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * New integrator - uses more history, but still no velocities
 *
 * This is effectively performing Richardson extrapolation on the standard Verlet!
 */
class VerletStock : public MultistepIntegrator {
public:
  VerletStock (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(4, _gravSys, _level, _dt)
  {}

  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    // note precision problems that will arise when positions are subtracted!
    // maybe the velocity version will be more precise?
    // just adding the parentheses increases accuracy!
    s[0].pos = s[3].pos + ((s[1].pos - s[4].pos)
             + 0.25*_dt*_dt* ( 5.0*s[1].acc
                              +2.0*s[2].acc
                              +5.0*s[3].acc));

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
  Hamming416 (NBodyGrav& _gravSys, const int _level, const double _dt) :
    MultistepIntegrator(4, _gravSys, _level, _dt)
  {}
  
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    s[0].pos = 2.0*s[2].pos - s[4].pos
             + (4.0*_dt*_dt/3.0) * ( s[1].acc
                                    +s[2].acc
                                    +s[3].acc);

    // get rid of oldest
    s.pop_back();
  }
};


/*
 * Hamming - pg 418
 *
 * A predictor using accelerations, velocities, and positions
 */
class Hamming418 {
public:
  Hamming418 (NBodyGrav& _gravSys, const int _level, const double _dt) :
    g(_gravSys),
    s(2, DynamicState(_level, 0))
  {
    s[0].step = 0;
    s[1].step = -1;
    // set initial conditions
    s[0].pos = g.initPosit();
    s[0].vel = g.initVeloc();
    // and set the previous state
    RK4 r(g,0);
    for (int istep=1; istep<2; ++istep) {
      for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
      s[istep].pos = r.getPosition();
      s[istep].vel = r.getVelocity();
      s[istep].acc = g.getAccel(s[istep].pos);
    }
  }
  
  ~Hamming418 () {};

  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration
    // first, use Adams Moulton to find the new velocity
    // NEED SOMETHING BETTER HERE!
    //s[1].vel = s[2].vel
    //         + (_dt/2.0)      * (     s[1].acc +     s[2].acc);
    // then use the Hamming method to find the position
    s[0].pos = s[1].pos
             + (_dt/2.0)      * (    -s[1].vel + 3.0*s[2].vel)
             + (_dt*_dt/12.0) * (17.0*s[1].acc + 7.0*s[2].acc);
    // this is horrible!
    s[0].vel = (384.0/_dt)    * (     s[2].pos - s[1].pos)
             + (1.0)          * (312.*s[1].vel + 73.*s[2].vel)
             + (-1.0*_dt)     * (110.*s[1].acc + 8.0*s[2].acc);

    // get rid of oldest
    s.pop_back();
  }

  ArrayXd getPosition () {
    return s[0].pos;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  vector<DynamicState> s;
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

  // find the "exact" solution for AB4 - use this as the exact solution for everyone
  double time = 0.0;
  int maxSteps = 100000;
  double dt = 10.0 / maxSteps;
  AB5 exact(s,0,dt);
  cout << "'Exact' solution is from running " << maxSteps << " steps of AB5 at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    time += dt;
  }

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
    VerletStock ves(s,0,dt);
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
      ves.stepForward(dt);
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

    ArrayXd sSolution = ves.getPosition();
    //cout << "  VerletStock solution: " << sSolution.segment(0,4).transpose() << endl;
    //cout << "  Error in Verlet-Stock is " << exact.getError(sSolution) << endl;

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
