#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
 * Compile with
 *
 * g++ -o verlet2 -std=c++11 -I/usr/include/eigen3 verlet2.cpp
 *
 * Copyright 2016 Mark J. Stock, markjstock@gmail.com
 */


/*
 * A pure sine wave, period=2pi, phase=0
 */
class SineWave {
public:
  SineWave(const int _num) {
    num = _num;
    //period.resize(num);
    //phase.resize(num);
    period = 2.5 + 1.5*ArrayXd::Random(num);
    phase = 3.1415927 * ArrayXd::Random(num);
    //cout << "Periods " << period.transpose() << endl;
    //for (int i=0; i<num; ++i) {
    //  period[i] = 1.0 + 3.0 * i / static_cast<double>(num);
    //  phase[i] = 2.0 * 3.1415927 * i / static_cast<double>(num);
    //}
  };

  ~SineWave() {};

  ArrayXd value(double time) {
    ArrayXd newVal;
    newVal.resize(num);
    for (int i=0; i<num; ++i) {
      newVal[i] = sin(time/period[i] + phase[i]);
    }
    return newVal;
  }

  ArrayXd veloc(double time) {
    ArrayXd newVal;
    newVal.resize(num);
    for (int i=0; i<num; ++i) {
      newVal[i] = cos(time/period[i] + phase[i]) / period[i];
    }
    return newVal;
  }

  ArrayXd accel(double time) {
    ArrayXd newVal;
    newVal.resize(num);
    for (int i=0; i<num; ++i) {
      newVal[i] = -sin(time/period[i] + phase[i]) / pow(period[i],2);
    }
    return newVal;
  }

private:
  int num;
  ArrayXd period;
  ArrayXd phase;
};


/*
 * Your basic 1-step Euler integrator
 */
class Euler {
public:
  Euler (ArrayXd _p, ArrayXd _v) :
    pos(_p),
    vel(_v)
  {
    // acceleration stays unassigned
    acc.resize(pos.size());
  }
  
  ~Euler () {};

  void stepForward (double _dt, ArrayXd _newAcc) {
    acc = _newAcc;
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
  ArrayXd pos;
  ArrayXd vel;
  ArrayXd acc;
};


/*
 * Standard Verlet (non-velocity) integrator
 */
class Verlet {
public:
  Verlet (ArrayXd _pnow, ArrayXd _pprev) {
    // set two positions
    pos[0] = _pprev;
    pos[1] = _pnow;
  }
  
  ~Verlet () {};

  void stepForward (double _dt, ArrayXd _newAcc) {
    ArrayXd newPos = 2.0*pos[1] - pos[0] + _dt*_dt*_newAcc;
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
  ArrayXd pos[2];
};


/*
 * New integrator - uses more history, but still no velocities
 */
class VerletStock {
public:
  VerletStock (ArrayXd _p0, ArrayXd _pm1, ArrayXd _pm2, ArrayXd _pm3,
               ArrayXd _am1, ArrayXd _am2) {
    // set two positions
    pos[0] = _pm3;
    pos[1] = _pm2;
    pos[2] = _pm1;
    pos[3] = _p0;
    acc[0] = _am2;
    acc[1] = _am1;
  }
  
  ~VerletStock () {};

  void stepForward (double _dt, ArrayXd _newAcc) {
    // compute the new position
    ArrayXd newPos = pos[3]
                  + pos[1]
                  - pos[0]
                  + 0.25*_dt*_dt* ( 5.0*_newAcc
                                   +2.0*acc[1]
                                   +5.0*acc[0]);
    // shift all values down the array
    pos[0] = pos[1];
    pos[1] = pos[2];
    pos[2] = pos[3];
    pos[3] = newPos;
    acc[0] = acc[1];
    acc[1] = _newAcc;
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
  ArrayXd pos[4];
  ArrayXd acc[2];
};


int main () {

  // iterate a sine wave for a few steps
  SineWave s(1000);

  // iterate a gravitational n-body system for a few steps
  //NBodyGrav g(100);

  double time = 0.0;
  double dt = 1.0;
  int maxSteps = 1000;
  cout << "Running " << maxSteps << " steps at dt= " << dt << endl;

  // initialize integrators
  Euler e(s.value(time), s.veloc(time));
  //cout << time;
  //cout << " " << e.getPosition().transpose();

  Verlet ve(s.value(time), s.value(time-dt));
  //cout << " " << ve.getPosition();

  VerletStock ves(s.value(time), s.value(time-dt), s.value(time-2.0*dt),
                  s.value(time-3.0*dt), s.accel(time-dt), s.accel(time-2.0*dt));
  //cout << " " << ves.getPosition();
  //cout << endl;

  // integrate forward
  for (int i=0; i<maxSteps; ++i) {
    // find the new acceleration
    ArrayXd newAccel = s.accel(time);

    // take the forward step
    e.stepForward(dt, newAccel);
    ve.stepForward(dt, newAccel);
    ves.stepForward(dt, newAccel);
    time += dt;

    // write out the new position 
    //cout << time << " " << e.getPosition().transpose();
    //cout << " " << ve.getPosition();
    //cout << " " << ves.getPosition()
    //cout << endl;
  }

  // how does each perform?
  ArrayXd solution = s.value(time);
  //cout << "Theoretical solution is:" << endl;
  //cout << time << " " << solution.transpose() << endl;
  cout << "Error in Euler is " << e.getError(solution) << endl;
  cout << "Error in Verlet is " << ve.getError(solution) << endl;
  cout << "Error in Verlet-Stock is " << ves.getError(solution) << endl;

  exit(0);
}
