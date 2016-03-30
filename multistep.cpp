#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
 * Compile with
 *
 * g++ -o verlet4 -std=c++11 -I/usr/include/eigen3 verlet4.cpp
 *
 * Copyright 2016 Mark J. Stock, markjstock@gmail.com
 *
 * Compute "true" solution once, using RK4, and compare all results to it
 *
 */


/*
 * Must make a virtual superclass for integrable systems
 */
//virtual class IntegrableSystem {
//class SineWave : IntegrableSystem {


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
  Euler (NBodyGrav& _gravSys) :
    g(_gravSys)
  {
    pos = g.initPosit();
    vel = g.initVeloc();
  }
  
  ~Euler () {};

  void stepForward (double _dt) {
    acc = g.getAccel(pos);
    pos = pos + _dt*vel + 0.5*_dt*_dt*acc;
    vel = vel + _dt*acc;
  }

  ArrayXd getPosition () {
    return pos;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-pos;
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  ArrayXd pos;
  ArrayXd vel;
  ArrayXd acc;
};


/*
 * Runge-Kutta 4th order
 */
class RK4 {
public:
  RK4 (NBodyGrav& _gravSys) :
    g(_gravSys)
  {
    pos = g.initPosit();
    vel = g.initVeloc();
  }
  
  ~RK4 () {};

  void stepForwardTake1 (double _dt) {
    double hdt = 0.5*_dt;
    // the k1, k2... are going to be acceleration estimates
    ArrayXd p1 = pos;
    ArrayXd k1 = g.getAccel(p1);

    ArrayXd p2 = p1 + hdt* (vel + 0.5*hdt*k1);
    ArrayXd k2 = g.getAccel(p2);

    ArrayXd p3 = p1 + hdt* (vel + 0.5*hdt*k2);
    ArrayXd k3 = g.getAccel(p3);

    ArrayXd p4 = p1 + _dt* (vel + 0.5*_dt*k3);
    ArrayXd k4 = g.getAccel(p4);

    // midpoint method uses beginning and ending velocities
    pos += 0.5 * _dt * vel;
    vel += _dt * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
    pos += 0.5 * _dt * vel;
  }

  void stepForwardTake2 (double _dt) {
    double hdt = 0.5*_dt;
    // the k1, k2... are going to be velocity estimates
    ArrayXd p1 = pos;
    ArrayXd a1 = g.getAccel(p1);
    ArrayXd v1 = vel + hdt*a1;

    ArrayXd p2 = p1 + hdt* (vel + 0.5*hdt*a1);
    ArrayXd a2 = g.getAccel(p2);
    ArrayXd v2 = vel + hdt*a2;

    ArrayXd p3 = p1 + hdt* (vel + 0.5*hdt*a2);
    ArrayXd a3 = g.getAccel(p3);
    ArrayXd v3 = vel + hdt*a3;

    ArrayXd p4 = p1 + _dt* (vel + 0.5*_dt*a3);
    ArrayXd a4 = g.getAccel(p4);
    ArrayXd v4 = vel + _dt*a4;

    // find mean of velocities
    pos += 0.5 * _dt * vel;
    vel = (v1 + 2.0*v2 + 2.0*v3 + v4) / 6.0;
    pos += 0.5 * _dt * vel;
  }

  void stepForward (double _dt) {
    double hdt = 0.5*_dt;
    // the k1, k2... are going to be acceleration estimates
    ArrayXd p1 = pos;
    ArrayXd k1 = g.getAccel(p1);

    ArrayXd p2 = p1 + hdt* (vel + 0.5*hdt*k1);
    ArrayXd k2 = g.getAccel(p2);

    ArrayXd p3 = p1 + hdt* (vel + 0.5*hdt*k2);
    ArrayXd k3 = g.getAccel(p3);

    ArrayXd p4 = p1 + _dt* (vel + 0.5*_dt*k3);
    ArrayXd k4 = g.getAccel(p4);

    // midpoint method uses beginning and ending velocities
    acc = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
    pos += _dt * (vel + 0.5*_dt*acc);
    vel += _dt * acc;
  }

  ArrayXd getPosition () {
    return pos;
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-pos;
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  ArrayXd pos;
  ArrayXd vel;
  ArrayXd acc;
};


/*
 * Standard Verlet (non-velocity) integrator
 */
class Verlet {
public:
  Verlet (NBodyGrav& _gravSys, double _dt) :
    g(_gravSys)
  {
    // set two positions
    pos[1] = g.initPosit();

    // need to find the previous position!
    RK4 rk(g);
    for (int i=0; i<100; ++i) {
      rk.stepForward(-0.01 * _dt);
    }
    pos[0] = rk.getPosition();
  }
  
  ~Verlet () {};

  void stepForward (double _dt) {
    ArrayXd newAcc = g.getAccel(pos[1]);
    ArrayXd newPos = 2.0*pos[1] - pos[0] + _dt*_dt*newAcc;
    pos[0] = pos[1];
    pos[1] = newPos;
  }

  ArrayXd getPosition () {
    return pos[1];
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  ArrayXd pos[2];
};


/*
 * New integrator - uses more history, but still no velocities
 */
class VerletStock {
public:
  VerletStock (NBodyGrav& _gravSys, double _dt) :
    g(_gravSys)
  {
    // set current position and acceleration
    pos[3] = g.initPosit();

    // need to find the previous position and accel
    RK4 rk(g);
    for (int i=0; i<100; ++i) {
      rk.stepForward(-0.01 * _dt);
    }
    //rk.stepForward(-1.0 * _dt);
    pos[2] = rk.getPosition();
    acc[1] = g.getAccel(pos[2]);

    // and two positions previous
    for (int i=0; i<100; ++i) {
      rk.stepForward(-0.01 * _dt);
    }
    //rk.stepForward(-1.0 * _dt);
    pos[1] = rk.getPosition();
    acc[0] = g.getAccel(pos[1]);

    for (int i=0; i<100; ++i) {
      rk.stepForward(-0.01 * _dt);
    }
    //rk.stepForward(-1.0 * _dt);
    pos[0] = rk.getPosition();
  }
  
  ~VerletStock () {};

  void stepForward (double _dt) {
    ArrayXd newAcc = g.getAccel(pos[3]);

    // compute the new position
    ArrayXd newPos = pos[3]
                   + pos[1]
                   - pos[0]
                   + 0.25*_dt*_dt* ( 5.0*newAcc
                                    +2.0*acc[1]
                                    +5.0*acc[0]);
    // shift all values down the array
    pos[0] = pos[1];
    pos[1] = pos[2];
    pos[2] = pos[3];
    pos[3] = newPos;
    acc[0] = acc[1];
    acc[1] = newAcc;
  }

  ArrayXd getPosition () {
    return pos[3];
  }

  double getError (ArrayXd _trueSolution) {
    ArrayXd temp = _trueSolution-getPosition();
    double normsq = temp.matrix().norm() / sqrt(temp.size());
    return(normsq);
  }

private:
  NBodyGrav& g;
  ArrayXd pos[4];
  ArrayXd acc[2];
};


// Perform Richardson Extrapolation?

// Create a system and an integrator
int main () {

  // iterate a sine wave for a few steps
  //SineWave s(1000);

  // iterate a gravitational n-body system for a few steps
  NBodyGrav s(100);

  double time = 0.0;
  int maxSteps = 1000;
  double dt = 10.0 / maxSteps;
  cout << "Running " << maxSteps << " steps at dt= " << dt << endl;

  // initialize integrators
  Euler e(s);
  RK4 r(s);
  //cout << time;
  //cout << " " << e.getPosition().transpose();

  Verlet ve(s, dt);
  //cout << " " << ve.getPosition();

  VerletStock ves(s, dt);
  //VerletStock ves(s.value(time), s.value(time-dt), s.value(time-2.0*dt),
  //                s.value(time-3.0*dt), s.accel(time-dt), s.accel(time-2.0*dt));
  //cout << " " << ves.getPosition();
  //cout << endl;

  // integrate forward
  for (int i=0; i<maxSteps; ++i) {

    // take the forward step
    e.stepForward(dt);
    r.stepForward(dt);
    ve.stepForward(dt);
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
  ArrayXd rSolution = r.getPosition();
  cout << "RK4 solution: " << rSolution.segment(0,4).transpose() << endl;
  ArrayXd vSolution = ve.getPosition();
  cout << "Verlet solution: " << vSolution.segment(0,4).transpose() << endl;
  ArrayXd sSolution = ves.getPosition();
  cout << "VerletStock solution: " << sSolution.segment(0,4).transpose() << endl;

  // find the "exact" solution?
  maxSteps = 10000;
  dt = 10.0 / maxSteps;

  // find the "exact" solution for RK4 - use this as the exact solution for everyone
  time = 0.0;
  RK4 ra(s);
  cout << "'Exact' solution is from running " << maxSteps << " steps of RK4 at dt= " << dt << endl;
  for (int i=0; i<maxSteps; ++i) {
    ra.stepForward(dt);
    time += dt;
  }

  cout << "Error in Euler is " << ra.getError(eSolution) << endl;
  cout << "Error in RK4 is " << ra.getError(rSolution) << endl;
  cout << "Error in Verlet is " << ra.getError(vSolution) << endl;
  cout << "Error in Verlet-Stock is " << ra.getError(sSolution) << endl;

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
