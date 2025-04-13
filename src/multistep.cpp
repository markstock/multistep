/*
 * multistep - a program to test forward advection schemes
 *
 * Copyright 2016,22,25 Mark J. Stock, markjstock@gmail.com
 */

#include "DynamicState.hpp"

#include "SineWave.hpp"
#include "SpringMass.hpp"
#include "LennardJones.hpp"
#include "NBodyGrav3D.hpp"
#include "NBodyVort2D.hpp"
#include "Projectiles3D.hpp"
//#include "NBodyVort3D.hpp"

#include "ForwardIntegrator.hpp"
#include "MultistageIntegrator.hpp"
#include "MultistepIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


// use this for sine waves
//#define TEMPLATEVAR double

// use this for everything else
#define TEMPLATEVAR Eigen::ArrayXd
// this is the same as Eigen::Array< double, Dynamic, 1 >.

// can we improve precision with this?
// #include <quadmath.h>
// Eigen::Array< __float128, Dynamic, 1 >.

// Create a system and an integrator
int main (int argc, char *argv[]) {

  constexpr int32_t numInteg = 16;
  std::string integ[numInteg] = {
    "Euler",	// 0  Euler
    "RK2-R",	// 1  RK2 Ralston
    "RK2-H",	// 2  RK2 Heun
    "AB2",		// 3  AB2
    "Verlet",	// 4  Verlet
    "RK3-R",	// 5  RK3 Ralston
    "RK3-H",	// 6  RK3 Heun
    "RK3-K",	// 7  RK3 Kutta
    "AB3",		// 8  AB3+AM3
    "RK4",		// 9  RK4 Classic
    "RK4-TE",	// 10 RK4 3/8ths rule
    "AB4",		// 11 AB4+AM4
    "RichVer",	// 12 Richardson-Verlet
    "Ham416",	// 13 Hamming416
    "Ham418",	// 14 Hamming418
    "AB5"		// AB5+AM5
  };
  // use numStages to scale the time steps when equivalentWork is true
  const int32_t numStages[numInteg] = {1, 2, 2, 1, 1, 3, 3, 3, 1, 4, 4, 1, 1, 1, 1, 1};

  bool equivalentWork = false;
  bool dumpevery = false;
  bool doRK = false;
  bool doAB = false;
  bool doHam = false;
  bool doVer = false;
  bool doInteg[numInteg] = {false};
  int32_t onceMaxSteps = 60;
  bool onceThrough = false;

  // parse command-line arguments
  // input file can appear anywhere
  // processing options are performed in the order they appear
  //(void) strcpy(progname,argv[0]);
  //if (argc < 2) (void) Usage(progname,0);
  for (int i=1; i<argc; i++) {
    if (strncmp(argv[i], "-",1) == 0) {
      if (strncmp(argv[i], "-all", 4) == 0) {
        doInteg[0] = true;
        doRK = true;
        doAB = true;
        doVer = true;
        doHam = true;
      } else if (strncmp(argv[i], "-rk2", 4) == 0) {
        doInteg[1] = true;
        doInteg[2] = true;
      } else if (strncmp(argv[i], "-rk3", 4) == 0) {
        doInteg[5] = true;
        doInteg[6] = true;
        doInteg[7] = true;
      } else if (strncmp(argv[i], "-rk4", 4) == 0) {
        doInteg[9] = true;
        doInteg[10] = true;
      } else if (strncmp(argv[i], "-rk", 3) == 0) {
        doRK = true;
      } else if (strncmp(argv[i], "-ab", 4) == 0) {
        doAB = true;
      } else if (strncmp(argv[i], "-ver", 4) == 0) {
        doVer = true;
      } else if (strncmp(argv[i], "-ham", 4) == 0) {
        doHam = true;
      } else if (strncmp(argv[i], "-eul", 4) == 0) {
        doInteg[0] = true;
      } else if (strncmp(argv[i], "-steps", 2) == 0) {
        onceMaxSteps = atoi(argv[++i]);
        onceThrough = true;
      } else if (strncmp(argv[i], "-equiv", 3) == 0) {
        equivalentWork = true;
      } else if (strncmp(argv[i], "-noequiv", 4) == 0) {
        equivalentWork = false;
      } else {
        // nothing?
      }
    }
  }

  // define the dynamical system
  //VelocitySine s(10.0/9.25);
  //AccelerationSine s(10.0/9.25);
  //AccelerationSine s(10.0/1.25);
  //AccelerationSine s(2.0*M_PI);
  //AccelerationSine s(40.0);
  SpringMass s(100,10.0/9.25);
  //LennardJones s(100,1.0,0.02);
  //NBodyGrav3D s(5);
  //NBodyGrav3D s(7);
  //NBodyVort2D s(32);
  //Projectiles3D s(64);

  // adjust available integrators based on dynamical system
  if (s.hasDamping() and doVer) {
    std::cout << "System has damping, turning off Verlet integrators\n";
    doVer = false;
  }
  if (doVer and not s.hasAccel()) {
    std::cout << "System does not use acceleration, turning off Verlet integrators\n";
    doVer = false;
  }
  if (doHam and not s.hasAccel()) {
    std::cout << "System does not use acceleration, turning off Hamming integrators\n";
    doHam = false;
  }

  if (doRK) {
    doInteg[0] = true;
    doInteg[1] = true;
    doInteg[2] = true;
    doInteg[5] = true;
    doInteg[6] = true;
    doInteg[7] = true;
    doInteg[9] = true;
    doInteg[10] = true;
  }
  if (doAB) {
    doInteg[0] = true;
    doInteg[3] = true;
    doInteg[8] = true;
    doInteg[11] = true;
    doInteg[15] = true;
  }
  if (doVer) {
    doInteg[4] = true;
    doInteg[12] = true;
  }
  if (doHam) {
    doInteg[13] = true;
    doInteg[14] = true;
  }


  // each system has its own end time
  double endtime = s.getEndTime();

  // find the "exact" solution and save it for reuse
  TEMPLATEVAR exact = s.getExact(endtime);

  std::cout << "steps\t";
  for (int32_t i=0; i<numInteg; ++i) {
    if (doInteg[i]) std::cout << integ[i] << "\t\t";
  }
  std::cout << std::endl;

  // set precision for error output
  std::cout << std::setprecision(8);

  // integrate using the various methods (always divisible by 12, though!)
  // ensure that maxsteps is divisible by 12 (to accommodate RK 2,3,4)
  //for (int32_t maxSteps = 60; maxSteps < 2000000; maxSteps = 12*(int)(maxSteps*1.3/12.0)) {
  for (int32_t maxSteps = 60; maxSteps < 130000; maxSteps = 12*(int)(maxSteps*1.3/12.0)) {

    if (onceThrough) maxSteps = onceMaxSteps;

    const double dt = endtime / maxSteps;
    if (dumpevery) std::cout << "Running " << maxSteps << " steps at dt= " << dt << std::endl;

    std::cout << maxSteps;

    // initialize integrators
    if (doInteg[0]) {
      Euler<TEMPLATEVAR> e(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[0] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[0]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        e.stepForward(thisDt);
        if (dumpevery) std::cout << "e at t " << e.getTime() << " is " << e.getPosition() << "\n";
      }
      //TEMPLATEVAR eSolution = e.getPosition();
      //std::cout << "  Euler solution: " << eSolution << std::endl;
      //std::cout << "  Euler solution: " << eSolution << std::endl;
      //std::cout << "  Error in Euler is " << exact.getError(eSolution) << std::endl;
      std::cout << "\t" << e.getError(exact);
    }

    if (doInteg[1]) {
      RK2Ralston<TEMPLATEVAR> r2(s,0);		// Ralston's is better for larger timesteps only
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[1] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[1]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r2.stepForward(thisDt);
        if (dumpevery) std::cout << "rk2r at t " << r2.getTime() << " is " << r2.getPosition() << "\n";
      }
      //TEMPLATEVAR r2Solution = r2.getPosition();
      //std::cout << "  RK2 solution: " << r2Solution << std::endl;
      //std::cout << "  RK2 solution: " << r2Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK2 is " << exact.getError(r2Solution) << std::endl;
      std::cout << "\t" << r2.getError(exact);
    }

    if (doInteg[2]) {
      RK2Heun<TEMPLATEVAR> r2(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[2] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[2]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r2.stepForward(thisDt);
        if (dumpevery) std::cout << "rk2h at t " << r2.getTime() << " is " << r2.getPosition() << "\n";
      }
      //TEMPLATEVAR r2Solution = r2.getPosition();
      //std::cout << "  RK2 solution: " << r2Solution << std::endl;
      //std::cout << "  RK2 solution: " << r2Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK2 is " << exact.getError(r2Solution) << std::endl;
      std::cout << "\t" << r2.getError(exact);
    }

    if (doInteg[3]) {
      AB2<TEMPLATEVAR> a2(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[3] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[3]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        a2.stepForward(thisDt);
        if (dumpevery) std::cout << "ab2 at t " << a2.getTime() << " is " << a2.getPosition() << "\n";
      }
      //TEMPLATEVAR a2Solution = a2.getPosition();
      //std::cout << "  ABM2 solution: " << a2Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in ABM2 is " << exact.getError(a2Solution) << std::endl;
      std::cout << "\t" << a2.getError(exact);
    }

    if (doInteg[4]) {
      Verlet<TEMPLATEVAR> ve(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[4] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[4]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        ve.stepForward(thisDt);
        if (dumpevery) std::cout << "ve at t " << ve.getTime() << " is " << ve.getPosition() << "\n";
      }
      //TEMPLATEVAR vSolution = ve.getPosition();
      //std::cout << "  Verlet solution: " << vSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in Verlet is " << exact.getError(vSolution) << std::endl;
      std::cout << "\t" << ve.getError(exact);
    }

    if (doInteg[5]) {
      RK3Ralston<TEMPLATEVAR> r3(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[5] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[5]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r3.stepForward(thisDt);
        if (dumpevery) std::cout << "rk3 at t " << r3.getTime() << " is " << r3.getPosition() << "\n";
      }
      //TEMPLATEVAR r3Solution = r3.getPosition();
      //std::cout << "  RK3 solution: " << r3Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK3 is " << exact.getError(r3Solution) << std::endl;
      std::cout << "\t" << r3.getError(exact);
    }

    if (doInteg[6]) {
      RK3Heun<TEMPLATEVAR> r3(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[6] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[6]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r3.stepForward(thisDt);
        if (dumpevery) std::cout << "rk3 at t " << r3.getTime() << " is " << r3.getPosition() << "\n";
      }
      //TEMPLATEVAR r3Solution = r3.getPosition();
      //std::cout << "  RK3 solution: " << r3Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK3 is " << exact.getError(r3Solution) << std::endl;
      std::cout << "\t" << r3.getError(exact);
    }

    if (doInteg[7]) {
      RK3Kutta<TEMPLATEVAR> r3(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[7] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[7]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r3.stepForward(thisDt);
        if (dumpevery) std::cout << "rk3 at t " << r3.getTime() << " is " << r3.getPosition() << "\n";
      }
      //TEMPLATEVAR r3Solution = r3.getPosition();
      //std::cout << "  RK3 solution: " << r3Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK3 is " << exact.getError(r3Solution) << std::endl;
      std::cout << "\t" << r3.getError(exact);
    }

    if (doInteg[8]) {
      AB3<TEMPLATEVAR> a3(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[8] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[8]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        a3.stepForward(thisDt);
        if (dumpevery) std::cout << "ab3 at t " << a3.getTime() << " is " << a3.getPosition() << "\n";
      }
      //TEMPLATEVAR a3Solution = a3.getPosition();
      //std::cout << "  ABM3 solution: " << a3Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in ABM3 is " << exact.getError(a3Solution) << std::endl;
      std::cout << "\t" << a3.getError(exact);
    }

    if (doInteg[9]) {
      RK4<TEMPLATEVAR> r4(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[9] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[9]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r4.stepForward(thisDt);
        if (dumpevery) std::cout << "rk4 at t " << r4.getTime() << " is " << r4.getPosition() << "\n";
      }
      //TEMPLATEVAR r4Solution = r4.getPosition();
      //std::cout << "  RK4 solution: " << r4Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK4 is " << exact.getError(r4Solution) << std::endl;
      std::cout << "\t" << r4.getError(exact);
    }

    if (doInteg[10]) {
      RK4ter<TEMPLATEVAR> r4(s,0);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[10] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[10]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        r4.stepForward(thisDt);
        if (dumpevery) std::cout << "rk4 at t " << r4.getTime() << " is " << r4.getPosition() << "\n";
      }
      //TEMPLATEVAR r4Solution = r4.getPosition();
      //std::cout << "  RK4 solution: " << r4Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK4 is " << exact.getError(r4Solution) << std::endl;
      std::cout << "\t" << r4.getError(exact);
    }

    if (doInteg[11]) {
      AB4<TEMPLATEVAR> ab(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[11] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[11]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        ab.stepForward(thisDt);
        if (dumpevery) std::cout << "ab4 at t " << ab.getTime() << " is " << ab.getPosition() << "\n";
      }
      //TEMPLATEVAR abSolution = ab.getPosition();
      //std::cout << "  ABM4 solution: " << abSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in ABM4 is " << exact.getError(abSolution) << std::endl;
      std::cout << "\t" << ab.getError(exact);
    }

    if (doInteg[12]) {
      RichardsonVerlet<TEMPLATEVAR> rv(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[12] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[12]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        rv.stepForward(thisDt);
        if (dumpevery) std::cout << "rv at t " << rv.getTime() << " is " << rv.getPosition() << "\n";
      }
      //TEMPLATEVAR rvSolution = rv.getPosition();
      //std::cout << "  RichardsonVerlet solution: " << sSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RichardsonVerlet is " << exact.getError(sSolution) << std::endl;
      std::cout << "\t" << rv.getError(exact);
    }

    if (doInteg[13]) {
      Hamming416<TEMPLATEVAR> h6(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[13] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[13]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        h6.stepForward(thisDt);
      }
      //templatevar h6solution = h6.getposition();
      std::cout << "\t" << h6.getError(exact);
    }

    if (doInteg[14]) {
      Hamming418<TEMPLATEVAR> h8(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[14] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[14]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        h8.stepForward(thisDt);
      }
      //templatevar h8solution = h8.getposition();
      std::cout << "\t" << h8.getError(exact);
    }

    if (doInteg[15]) {
      AB5<TEMPLATEVAR> abb(s,0,dt);
      const int32_t thisSteps = equivalentWork ? maxSteps/numStages[15] : maxSteps;
      const double  thisDt    = equivalentWork ? numStages[15]*dt       : dt;
      for (int32_t i=0; i<thisSteps; ++i) {
        abb.stepForward(thisDt);
        if (dumpevery) std::cout << "ab5 at t " << abb.getTime() << " is " << abb.getPosition() << "\n";
      }
      //templatevar abbsolution = abb.getposition();
      std::cout << "\t" << abb.getError(exact);
    }

    std::cout << std::endl;
    if (onceThrough) break;
  }

  exit(0);
}
