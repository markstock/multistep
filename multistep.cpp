#include <iostream>
#include <math.h>

using namespace std;

/*
 * A pure sine wave, period=2pi, phase=0
 */
class SineWave {
public:
  SineWave() {};
  ~SineWave() {};

  double value(double time) {
    return sin(time);
  }
  double veloc(double time) {
    return cos(time);
  }
  double accel(double time) {
    return -sin(time);
  }

private:
};

/*
 * Your basic 1-step Euler integrator
 */
class Euler {
public:
  Euler (double _p, double _v) :
    pos(_p),
    vel(_v)
  {
    // acceleration stays unassigned
    acc = 0.0;
  }
  
  ~Euler () {};

  void stepForward (double _dt, double _newAcc) {
    acc = _newAcc;
    pos = pos + _dt*vel + 0.5*_dt*_dt*acc;
    vel = vel + _dt*acc;
  }

  double getPosition () {
    return pos;
  }

private:
  double pos;
  double vel;
  double acc;
};

/*
 * Standard Verlet (non-velocity) integrator
 */
class Verlet {
public:
  Verlet (double _pnow, double _pprev) {
    // set two positions
    pos[0] = _pprev;
    pos[1] = _pnow;
  }
  
  ~Verlet () {};

  void stepForward (double _dt, double _newAcc) {
    double newPos = 2.0*pos[1] - pos[0] + _dt*_dt*_newAcc;
    pos[0] = pos[1];
    pos[1] = newPos;
  }

  double getPosition () {
    return pos[1];
  }

private:
  double pos[2];
};


/*
 * New integrator - uses more history, but still no velocities
 */
class VerletStock {
public:
  VerletStock (double _p0, double _pm1, double _pm2, double _pm3,
               double _am1, double _am2) {
    // set two positions
    pos[0] = _pm3;
    pos[1] = _pm2;
    pos[2] = _pm1;
    pos[3] = _p0;
    acc[0] = _am2;
    acc[1] = _am1;
  }
  
  ~VerletStock () {};

  void stepForward (double _dt, double _newAcc) {
    // compute the new position
    double newPos = pos[3]
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

  double getPosition () {
    return pos[3];
  }

private:
  double pos[4];
  double acc[2];
};


int main () {

  // iterate a sine wave for a few steps
  SineWave s;
  double time = 0.0;
  double dt = 0.0001;
  int maxSteps = 1000;
  cout << time;


  // initialize integrators
  Euler e(s.value(time), s.veloc(time));
  cout << " " << e.getPosition();

  Verlet ve(s.value(time), s.value(time-dt));
  cout << " " << ve.getPosition();

  VerletStock ves(s.value(time), s.value(time-dt), s.value(time-2.0*dt),
                  s.value(time-3.0*dt), s.accel(time-dt), s.accel(time-2.0*dt));
  cout << " " << ves.getPosition();
  cout << endl;


  // integrate forward
  for (int i=0; i<maxSteps; ++i) {
    // find the new acceleration
    double newAccel = s.accel(time);

    // take the forward step
    e.stepForward(dt, newAccel);
    ve.stepForward(dt, newAccel);
    ves.stepForward(dt, newAccel);
    time += dt;

    // write out the new position 
    cout << time << " " << e.getPosition();
    cout << " " << ve.getPosition();
    cout << " " << ves.getPosition() << endl;
  }

  // how does each perform?
  cout << "Theoretical solution at t=" << time << " is " << s.value(time) << endl;
  cout << "Error in Euler is " << (fabs(s.value(time)-e.getPosition())/s.value(time)) << endl;
  cout << "Error in Verlet is " << (fabs(s.value(time)-ve.getPosition())/s.value(time)) << endl;
  cout << "Error in VerletStock is " << fabs((s.value(time)-ves.getPosition())/s.value(time)) << endl;

  exit(0);
}
