/*
 * multistep - a program to test forward advection schemes
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#include "DynamicState.hpp"

#include "NBodyVort2D.hpp"
#include "NBodyGrav3D.hpp"
#include "SineWave.hpp"
#include "SpringMass.hpp"
#include "LennardJones.hpp"

#include "ForwardIntegrator.hpp"
#include "MultistageIntegrator.hpp"
#include "MultistepIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


// Perform Richardson Extrapolation by creating multiple temporal resolution levels automatically?
// or is that Bulirsch-Stoer?

/*
 * Control class for Richardson Extrapolation and globally adaptive time stepping
 */
class RichardsonEuler {
public:
  RichardsonEuler (DynamicalSystem<Eigen::ArrayXd>& _system) {

  }

private:
  //Euler
};

// use this for sine waves
//#define TEMPLATEVAR double

// use this for NBodyGrav3D and NBodyVort2D
#define TEMPLATEVAR Eigen::ArrayXd

// Create a system and an integrator
int main () {

  const bool dumpevery = false;
  const bool doHamming = false;
  const double endtime = 10.0;

  // define the dynamical system
  //VelocitySine s(endtime/9.25);
  //AccelerationSine s(endtime/9.25);
  //SpringMass s(100,endtime/9.25);
  LennardJones s(100,1.0,0.02);
  //AccelerationSine s(endtime/1.25);
  //AccelerationSine s(2.0*M_PI);
  //AccelerationSine s(endtime*4);
  //NBodyGrav3D s(10);
  //NBodyVort2D s(10);

  // find the "exact" solution and save it for reuse
  TEMPLATEVAR exact = s.getExact(endtime);

  std::cout << "steps\t";
  std::cout << "Euler\t\t";
  std::cout << "RK2\t\t";
  std::cout << "AB2\t\t";
  if (s.hasAccel()) std::cout << "Verlet\t\t";
  std::cout << "RK3\t\t";
  std::cout << "RK4\t\t";
  std::cout << "AB4\t\t";
  if (s.hasAccel()) std::cout << "Verlet4\t\t";
  if (s.hasAccel() and doHamming) std::cout << "Ham416\t\t";
  if (s.hasAccel() and doHamming) std::cout << "Ham418\t\t";
  std::cout << "AB5";
  std::cout << std::endl;

  // set precision for error output
  std::cout << std::setprecision(8);

  // integrate using the various methods
  for (int32_t maxSteps = 60; maxSteps < 130000; maxSteps = 12*(int)(maxSteps*1.3/12.0)) {
  //for (int32_t maxSteps = 60; maxSteps < 130000; maxSteps *= 2) {
  //for (int32_t maxSteps = 12; maxSteps < 15; maxSteps *= 2) {
  //for (int32_t maxSteps = 100; maxSteps < 105; maxSteps *= 2) {

    // ensure that maxsteps is divisible by 12 (to accomodate RK 2,3,4)

    const double dt = endtime / maxSteps;
    if (dumpevery) std::cout << "Running " << maxSteps << " steps at dt= " << dt << std::endl;

    std::cout << maxSteps;

    // initialize integrators
    {
      Euler<TEMPLATEVAR> e(s,0);
      for (int32_t i=0; i<maxSteps; ++i) {
        //std::cout << "\nvalue is " << e.getPosition();
        e.stepForward(dt);
        if (dumpevery) std::cout << "e at t " << e.getTime() << " is " << e.getPosition() << "\n";
      }
      //TEMPLATEVAR eSolution = e.getPosition();
      //std::cout << "  Euler solution: " << eSolution << std::endl;
      //std::cout << "  Euler solution: " << eSolution << std::endl;
      //std::cout << "  Error in Euler is " << exact.getError(eSolution) << std::endl;
      std::cout << "\t" << e.getError(exact);
    }

    {
      RK2Ralston<TEMPLATEVAR> r2(s,0);		// Ralston's is better for larger timesteps only
      for (int32_t i=0; i<maxSteps; ++i) {
        if (i%2 == 0) {
          r2.stepForward(2.0*dt);
          if (dumpevery) std::cout << "rk2 at t " << r2.getTime() << " is " << r2.getPosition() << "\n";
        }
      }
      //TEMPLATEVAR r2Solution = r2.getPosition();
      //std::cout << "  RK2 solution: " << r2Solution << std::endl;
      //std::cout << "  RK2 solution: " << r2Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK2 is " << exact.getError(r2Solution) << std::endl;
      std::cout << "\t" << r2.getError(exact);
    }

    {
      AB2<TEMPLATEVAR> a2(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        a2.stepForward(dt);
        if (dumpevery) std::cout << "ab2 at t " << a2.getTime() << " is " << a2.getPosition() << "\n";
      }
      //TEMPLATEVAR a2Solution = a2.getPosition();
      //std::cout << "  ABM2 solution: " << a2Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in ABM2 is " << exact.getError(a2Solution) << std::endl;
      std::cout << "\t" << a2.getError(exact);
    }

    if (s.hasAccel()) {
      Verlet<TEMPLATEVAR> ve(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        ve.stepForward(dt);
        if (dumpevery) std::cout << "ve at t " << ve.getTime() << " is " << ve.getPosition() << "\n";
      }
      //TEMPLATEVAR vSolution = ve.getPosition();
      //std::cout << "  Verlet solution: " << vSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in Verlet is " << exact.getError(vSolution) << std::endl;
      std::cout << "\t" << ve.getError(exact);
    }

    {
      RK3<TEMPLATEVAR> r3(s,0);
      for (int32_t i=0; i<maxSteps; ++i) {
        if (i%3 == 0) {
          r3.stepForward(3.0*dt);
          if (dumpevery) std::cout << "rk3 at t " << r3.getTime() << " is " << r3.getPosition() << "\n";
        }
      }
      //TEMPLATEVAR r3Solution = r3.getPosition();
      //std::cout << "  RK3 solution: " << r3Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK3 is " << exact.getError(r3Solution) << std::endl;
      std::cout << "\t" << r3.getError(exact);
    }

    {
      RK4<TEMPLATEVAR> r4(s,0);
      for (int32_t i=0; i<maxSteps; ++i) {
        if (i%4 == 0) {
          r4.stepForward(4.0*dt);
          if (dumpevery) std::cout << "rk4 at t " << r4.getTime() << " is " << r4.getPosition() << "\n";
        }
      }
      //TEMPLATEVAR r4Solution = r4.getPosition();
      //std::cout << "  RK4 solution: " << r4Solution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RK4 is " << exact.getError(r4Solution) << std::endl;
      std::cout << "\t" << r4.getError(exact);
    }

    {
      AB4<TEMPLATEVAR> ab(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        ab.stepForward(dt);
        if (dumpevery) std::cout << "ab4 at t " << ab.getTime() << " is " << ab.getPosition() << "\n";
      }
      //TEMPLATEVAR abSolution = ab.getPosition();
      //std::cout << "  ABM4 solution: " << abSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in ABM4 is " << exact.getError(abSolution) << std::endl;
      std::cout << "\t" << ab.getError(exact);
    }

    if (s.hasAccel()) {
      RichardsonVerlet<TEMPLATEVAR> rv(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        rv.stepForward(dt);
        if (dumpevery) std::cout << "rv at t " << rv.getTime() << " is " << rv.getPosition() << "\n";
      }
      //TEMPLATEVAR rvSolution = rv.getPosition();
      //std::cout << "  RichardsonVerlet solution: " << sSolution.segment(0,4).transpose() << std::endl;
      //std::cout << "  Error in RichardsonVerlet is " << exact.getError(sSolution) << std::endl;
      std::cout << "\t" << rv.getError(exact);
    }

    if (s.hasAccel() and doHamming) {
      Hamming416<TEMPLATEVAR> h6(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        h6.stepForward(dt);
      }
      //templatevar h6solution = h6.getposition();
      std::cout << "\t" << h6.getError(exact);
    }

    if (s.hasAccel() and doHamming) {
      Hamming418<TEMPLATEVAR> h8(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        h8.stepForward(dt);
      }
      //templatevar h8solution = h8.getposition();
      std::cout << "\t" << h8.getError(exact);
    }

    {
      AB5<TEMPLATEVAR> abb(s,0,dt);
      for (int32_t i=0; i<maxSteps; ++i) {
        abb.stepForward(dt);
        if (dumpevery) std::cout << "ab5 at t " << abb.getTime() << " is " << abb.getPosition() << "\n";
      }
      //templatevar abbsolution = abb.getposition();
      std::cout << "\t" << abb.getError(exact);
    }

    std::cout << std::endl;
  }

  exit(0);
}
