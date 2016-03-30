#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
 * Compile with
 *
 * g++ -o verlet5 -Ofast -std=c++11 -I/usr/include/eigen3 verlet5.cpp
 *
 * Copyright 2016 Mark J. Stock, markjstock@gmail.com
 *
 * Attempt to use Richardson extrapolation to refine solution, so need state objects
 *
 */


/*
 * Must make a virtual superclass for integrable systems
 */
//virtual class IntegrableSystem {
//class SineWave : IntegrableSystem {


/*
 * A class to contain one state (one set of changing variables)
 */
class DynamicState {
public: DynamicState(int _level, int _step) :
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
};


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
  ArrayXd getAccel(ArrayXd pos) {
    //cout << "pos in getAccel is " << pos.transpose() << endl;
    ArrayXd newVal = ArrayXd::Zero(3*num);
    //cout << "newVal in getAccel is " << newVal.transpose() << endl;
    for (int i=0; i<num; ++i) {
      // new accelerations on particle i
      Vector3d temp = pos.segment(3*i,3);
      //cout << "  particle i is at " << temp.transpose() << endl;
      Vector3d newAcc(0.0, 0.0, 0.0);
      for (int j=0; j<num; ++j) {
        // 20 flops
        // the influence of particle j
        //cout << "  particle j is at " << pos.segment(3*i,3).transpose() << endl;
        Vector3d dx = pos.segment(3*j,3) - pos.segment(3*i,3);
        double invdist = 1.0/(dx.norm()+radiusSquared(j));
        newAcc += dx * (mass(j) * invdist * invdist * invdist);
        //cout << "  influence of " << j << " on " << i << " is " << newAcc.transpose() << endl;
      }
      newVal.segment(3*i,3) = newAcc;
    }
    return newVal;
  }

private:
  int num;
  ArrayXd mass;
  ArrayXd radiusSquared;
  ArrayXd initPos;
  ArrayXd initVel;
};


/*
 * Your basic 1-step Euler integrator
 */
class Euler {
public:
  Euler (NBodyGrav& _gravSys, int _level) :
    g(_gravSys),
    curr(_level, 0),
    last(_level,-1)
  {
    // set initial conditions
    curr.pos = g.initPosit();
    curr.vel = g.initVeloc();
    // now curr has pos and vel but not acc
  }
  
  ~Euler () {};

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    // find the new force
    curr.acc = g.getAccel(curr.pos);
    // copy current state to past (deep copy?)
    last = curr;
    // perform forward integration
    curr.step++;
    curr.pos = last.pos + _dt*last.vel + 0.5*_dt*_dt*last.acc;
    curr.vel = last.vel + _dt*last.acc;
    // now curr has pos and vel but not acc
  }

  ArrayXd getPosition () {
    return curr.pos;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  DynamicState last;
  DynamicState curr;
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
 * Runge-Kutta 2nd order
 */
class RK2 {
public:
  RK2 (NBodyGrav& _gravSys, int _level) :
    g(_gravSys),
    curr(_level, 0),
    last(_level,-1)
  {
    // set initial conditions
    curr.pos = g.initPosit();
    curr.vel = g.initVeloc();
    // now curr has pos and vel but not acc
  }
  
  ~RK2 () {};

  // Generalizable to Heun's or Ralson's Methods
  void stepForward (const double _dt) {
    // errors at 800 steps are  0.5:1.77703e-05  1.0:1.78663e-05  2/3:1.77534e-05
    const double alpha = 2.0/3.0;
    curr.acc = g.getAccel(curr.pos);

    // first step: find acceleration at initial state
    ArrayXd pos1 = curr.pos;
    ArrayXd acc1 = curr.acc;
    ArrayXd vel1 = curr.vel;

    // second step: project forward a half step using that acceleration
    const double hdt = alpha*_dt;
    ArrayXd pos2 = pos1 + hdt*vel1;
    ArrayXd acc2 = g.getAccel(pos2);
    ArrayXd vel2 = vel1 + hdt*acc1;

    // copy current state to past
    last = curr;
    // perform forward integration
    curr.step++;
    // position updates via weighted average velocity
    const double oo2a = 1.0 / (2.0*alpha);
    curr.pos += _dt * ((1.0-oo2a)*vel1 + oo2a*vel2);
    // velocity updates via weighted average acceleration
    curr.vel += _dt * ((1.0-oo2a)*acc1 + oo2a*acc2);
  }

  ArrayXd getPosition () {
    return curr.pos;
  }

  ArrayXd getVelocity () {
    return curr.vel;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  DynamicState last;
  DynamicState curr;
};


/*
 * Runge-Kutta 4th order
 */
class RK4 {
public:
  RK4 (NBodyGrav& _gravSys, int _level) :
    g(_gravSys),
    curr(_level, 0),
    last(_level,-1)
  {
    // set initial conditions
    curr.pos = g.initPosit();
    curr.vel = g.initVeloc();
    // now curr has pos and vel but not acc
  }
  
  ~RK4 () {};

  // Most write-ups of this are incorrect! Does nobody edit their shit?
  // Note that this could be improved using the 3/8 rule, see wikipedia
  void stepForward (const double _dt) {
    curr.acc = g.getAccel(curr.pos);

    // first step: find acceleration at initial state
    ArrayXd pos1 = curr.pos;
    ArrayXd acc1 = curr.acc;
    ArrayXd vel1 = curr.vel;

    // second step: project forward a half step using that acceleration
    const double hdt = 0.5*_dt;
    ArrayXd pos2 = pos1 + hdt*vel1;
    ArrayXd acc2 = g.getAccel(pos2);
    ArrayXd vel2 = vel1 + hdt*acc1;

    // third step: project forward a half step from initial using the new acceleration
    ArrayXd pos3 = pos1 + hdt*vel2;
    ArrayXd acc3 = g.getAccel(pos3);
    ArrayXd vel3 = vel1 + hdt*acc2;

    // fourth step: project forward a full step from initial using the newest acceleration
    ArrayXd pos4 = pos1 + _dt*vel3;
    ArrayXd acc4 = g.getAccel(pos4);
    ArrayXd vel4 = vel1 + _dt*acc3;

    // copy current state to past (deep copy?)
    last = curr;
    // perform forward integration
    curr.step++;
    // position updates via weighted average velocity
    curr.pos += _dt * (vel1 + 2.0*vel2 + 2.0*vel3 + vel4) / 6.0;
    // velocity updates via weighted average acceleration
    curr.vel += _dt * (acc1 + 2.0*acc2 + 2.0*acc3 + acc4) / 6.0;
  }

  ArrayXd getPosition () {
    return curr.pos;
  }

  ArrayXd getVelocity () {
    return curr.vel;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  DynamicState last;
  DynamicState curr;
};


/*
 * Adams-Bashforth plus Adams-Moulton
 */
class AB2 {
public:
  AB2 (NBodyGrav& _gravSys, int _level, double _dt) :
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
    r.stepForward(-1.0 * _dt);
    s[1].pos = r.getPosition();
    s[1].vel = r.getVelocity();
    s[1].acc = g.getAccel(s[1].pos);
  }
  
  ~AB2 () {};

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration: AB2 for velocity, AM2 for position
    s[0].vel = s[1].vel + 0.5 * _dt * (3.0*s[1].acc - s[2].acc);
    s[0].pos = s[1].pos + 0.5 * _dt * (s[0].vel + s[1].vel);

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
 * Adams-Bashforth plus Adams-Moulton, all 4th order
 */
class AB4 {
public:
  AB4 (NBodyGrav& _gravSys, int _level, double _dt) :
    g(_gravSys),
    s(4, DynamicState(_level, 0))
  {
    s[0].step = 0;
    s[1].step = -1;
    s[2].step = -2;
    s[3].step = -3;
    // set initial conditions
    s[0].pos = g.initPosit();
    s[0].vel = g.initVeloc();
    // and set the previous state
    RK4 r(g,0);
    for (int istep=1; istep<4; ++istep) {
      for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
      s[istep].pos = r.getPosition();
      s[istep].vel = r.getVelocity();
      s[istep].acc = g.getAccel(s[istep].pos);
    }
  }
  
  ~AB4 () {};

  // always takes a position and velocity and turns it into a new position and velocity
  void stepForward (const double _dt) {
    s[0].acc = g.getAccel(s[0].pos);

    // add a new one to the head
    DynamicState newHead(s[0].level, s[0].step++);
    s.insert(s.begin(), newHead);

    // perform forward integration: AB4 for velocity, AM4 for position
    s[0].vel = s[1].vel +  _dt * (55.0*s[1].acc - 59.0*s[2].acc + 37.0*s[3].acc - 9.0*s[4].acc) / 24.0;
    s[0].pos = s[1].pos +  _dt * (9.0*s[0].vel + 19.0*s[1].vel - 5.0*s[2].vel + s[3].vel) / 24.0;

    // try AM5 for position - this is strangely worse
    //s[0].pos = s[1].pos +  _dt * (251.0*s[0].vel + 646.0*s[1].vel - 264.0*s[2].vel + 106.0*s[3].vel - 19.0*s[4].vel) / 720.0;

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
 * Standard Verlet (non-velocity) integrator
 */
class Verlet {
public:
  Verlet (NBodyGrav& _gravSys, int _level, double _dt) :
    g(_gravSys),
    s(2, DynamicState(_level, 0))
  {
    s[0].step = 0;
    s[1].step = -1;
    // set two positions
    s[0].pos = g.initPosit();

    // need to find the previous position, but not vel or acc
    RK4 r(g,0);
    for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
    s[1].pos = r.getPosition();
  }
  
  ~Verlet () {};

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
 * New integrator - uses more history, but still no velocities
 *
 * This is effectively performing Richardson extrapolation on the standard Verlet!
 */
class VerletStock {
public:
  VerletStock (NBodyGrav& _gravSys, int _level, double _dt) :
    g(_gravSys),
    s(4, DynamicState(_level, 0))
  {
    s[0].step = 0;
    s[1].step = -1;
    s[2].step = -2;
    s[3].step = -3;
    // set initial conditions
    s[0].pos = g.initPosit();
    // and set the previous state
    RK4 r(g,0);
    for (int istep=1; istep<4; ++istep) {
      for (int i=0; i<100; ++i) r.stepForward(-0.01 * _dt);
      s[istep].pos = r.getPosition();
      s[istep].acc = g.getAccel(s[istep].pos);
    }
  }
  
  ~VerletStock () {};

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

  // iterate a sine wave for a few steps
  //SineWave s(1000);

  // iterate a gravitational n-body system for a few steps
  NBodyGrav s(100);

  double time = 0.0;
  int maxSteps = 12800;
  double dt = 10.0 / maxSteps;
  cout << "Running " << maxSteps << " steps at dt= " << dt << endl;

  // initialize integrators
  Euler e(s,0);
  RK2 r2(s,0);
  AB2 a(s,0,dt);
  Verlet ve(s,0,dt);
  RK4 r(s,0);
  AB4 ab(s,0,dt);
  VerletStock ves(s,0,dt);
  //cout << time;
  //cout << " " << e.getPosition().transpose();

  //cout << " " << ve.getPosition();

  //VerletStock ves(s.value(time), s.value(time-dt), s.value(time-2.0*dt),
  //                s.value(time-3.0*dt), s.accel(time-dt), s.accel(time-2.0*dt));
  //cout << " " << ves.getPosition();
  //cout << endl;

  // integrate forward
  for (int i=0; i<maxSteps; ++i) {

    // take the forward step
    e.stepForward(dt);
    if (i%2 == 0) r2.stepForward(2.0*dt);
    a.stepForward(dt);
    ve.stepForward(dt);
    if (i%4 == 0) r.stepForward(4.0*dt);
    ab.stepForward(dt);
    ves.stepForward(dt);
    time += dt;

    // write out the new position 
    //cout << time << " " << e.getPosition().transpose();
    //cout << " " << ve.getPosition();
    //cout << " " << ves.getPosition()
    //cout << endl;
  }
  ArrayXd eSolution = e.getPosition();
  cout << "Euler solution: " << eSolution.segment(0,4).transpose() << endl;

  ArrayXd r2Solution = r2.getPosition();
  cout << "RK2 solution: " << r2Solution.segment(0,4).transpose() << endl;

  ArrayXd aSolution = a.getPosition();
  cout << "ABM2 solution: " << aSolution.segment(0,4).transpose() << endl;

  ArrayXd vSolution = ve.getPosition();
  cout << "Verlet solution: " << vSolution.segment(0,4).transpose() << endl;

  ArrayXd rSolution = r.getPosition();
  cout << "RK4 solution: " << rSolution.segment(0,4).transpose() << endl;

  ArrayXd abSolution = ab.getPosition();
  cout << "ABM4 solution: " << abSolution.segment(0,4).transpose() << endl;

  ArrayXd sSolution = ves.getPosition();
  cout << "VerletStock solution: " << sSolution.segment(0,4).transpose() << endl;

  // find the "exact" solution?
  maxSteps = 100000;
  dt = 10.0 / maxSteps;

  // find the "exact" solution for Euler - use this as the exact solution for everyone
/*
  time = 0.0;
  Euler exact(s,0);
  cout << "'Exact' solution is from running " << maxSteps << " steps of Euler at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    time += dt;
  }
  cout << "Error in Euler is " << exact.getError(eSolution) << endl;
*/

  // find the "exact" solution for RK4 - use this as the exact solution for everyone
/*
  time = 0.0;
  RK4 exact(s,0);
  cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    time += dt;
  }
*/

  // find the "exact" solution for AB4 - use this as the exact solution for everyone
  time = 0.0;
  AB4 exact(s,0,dt);
  cout << "'Exact' solution is from running " << maxSteps << " steps of AB4 at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    time += dt;
  }

  // find the "exact" solution for Verlet-Stock - use this as the exact solution for everyone
/*
  time = 0.0;
  VerletStock exact(s,0,dt);
  cout << "'Exact' solution is from running " << maxSteps << " steps of Verlet-Stock at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    time += dt;
  }
*/

  cout << "Error in Euler is " << exact.getError(eSolution) << endl;
  cout << "Error in RK2 is " << exact.getError(r2Solution) << endl;
  cout << "Error in ABM2 is " << exact.getError(aSolution) << endl;
  cout << "Error in Verlet is " << exact.getError(vSolution) << endl;
  cout << "Error in RK4 is " << exact.getError(rSolution) << endl;
  cout << "Error in ABM4 is " << exact.getError(abSolution) << endl;
  cout << "Error in Verlet-Stock is " << exact.getError(sSolution) << endl;

/*
  // for Euler
  time = 0.0;
  Euler ea(s);
  for (int i=0; i<maxSteps; ++i) {
    ea.stepForward(dt);
    time += dt;
  }
  cout << "Error in Euler is " << ea.getError(eSolution) << endl;

  // find the "exact" solution for Verlet
  time = 0.0;
  Verlet va(s, dt);
  for (int i=0; i<maxSteps; ++i) {
    va.stepForward(dt);
    time += dt;
  }
  cout << "Error in Verlet is " << va.getError(vSolution) << endl;

  // find the "exact" solution for Verlet-Stock
  time = 0.0;
  VerletStock vsa(s, dt);
  for (int i=0; i<maxSteps; ++i) {
    vsa.stepForward(dt);
    time += dt;
  }
  cout << "Error in Verlet-Stock is " << vsa.getError(sSolution) << endl;
*/

  exit(0);
}
