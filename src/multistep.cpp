/*
 * multistep - a program to test forward advection schemes
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#include "DynamicState.hpp"
#include "DynamicalSystem.hpp"
#include "NBodyVort2D.hpp"
#include "NBodyGrav3D.hpp"
#include "ForwardIntegrator.hpp"
#include "MultistageIntegrator.hpp"
#include "MultistepIntegrator.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <iostream>
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



// Create a system and an integrator
int main () {

  // iterate a gravitational n-body system for a few steps
  NBodyGrav3D s(100);
  //NBodyVort2D vort(100);

  // find the "exact" solution for AB4 - use this as the exact solution for everyone
  double time = 0.0;
  int32_t maxSteps = 10000;
  double dt = 10.0 / maxSteps;
  AB5<Eigen::ArrayXd> exact(s,0,dt);
  std::cout << "'Exact' solution is from running " << maxSteps << " steps of AB5 at dt= " << dt << std::endl;
  for (int32_t i=0; i<maxSteps; ++i) {
    exact.stepForward(dt);
    //ev.stepForward(dt);
    time += dt;
  }

  std::cout << "steps\t";
  std::cout << "Euler\t\t";
  std::cout << "RK2\t\t";
  std::cout << "AB2\t\t";
  std::cout << "Verlet\t\t";
  std::cout << "RK3\t\t";
  std::cout << "RK4\t\t";
  std::cout << "AB4\t\t";
  std::cout << "Verlet4\t\t";
  std::cout << "Ham416\t\t";
  //std::cout << "Ham418\t\t";
  std::cout << "AB5";
  std::cout << std::endl;

  // integrate using the various methods
  for (maxSteps = 12; maxSteps < 15000; maxSteps *= 2) {
    dt = 10.0 / maxSteps;
    //std::cout << "Running " << maxSteps << " steps at dt= " << dt << std::endl;
    time = 0.0;

    // initialize integrators
    Euler<Eigen::ArrayXd> e(s,0);
    RK2Ralston<Eigen::ArrayXd> r2(s,0);		// Ralston's is better for larger timesteps only
    AB2<Eigen::ArrayXd> a2(s,0,dt);
    Verlet<Eigen::ArrayXd> ve(s,0,dt);
    RK3<Eigen::ArrayXd> r3(s,0);
    RK4<Eigen::ArrayXd> r4(s,0);
    AB4<Eigen::ArrayXd> ab(s,0,dt);
    RichardsonVerlet<Eigen::ArrayXd> rv(s,0,dt);
    Hamming416<Eigen::ArrayXd> h6(s,0,dt);
    Hamming418<Eigen::ArrayXd> h8(s,0,dt);
    AB5<Eigen::ArrayXd> abb(s,0,dt);

    for (int32_t i=0; i<maxSteps; ++i) {

      // take the forward step
      e.stepForward(dt);
      if (i%2 == 0) r2.stepForward(2.0*dt);
      a2.stepForward(dt);
      ve.stepForward(dt);
      if (i%3 == 0) r3.stepForward(3.0*dt);
      if (i%4 == 0) r4.stepForward(4.0*dt);
      ab.stepForward(dt);
      rv.stepForward(dt);
      h6.stepForward(dt);
      h8.stepForward(dt);
      abb.stepForward(dt);
      time += dt;

    }

    Eigen::ArrayXd eSolution = e.getPosition();
    //std::cout << "  Euler solution: " << eSolution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in Euler is " << exact.getError(eSolution) << std::endl;

    Eigen::ArrayXd r2Solution = r2.getPosition();
    //std::cout << "  RK2 solution: " << r2Solution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in RK2 is " << exact.getError(r2Solution) << std::endl;

    Eigen::ArrayXd a2Solution = a2.getPosition();
    //std::cout << "  ABM2 solution: " << a2Solution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in ABM2 is " << exact.getError(a2Solution) << std::endl;

    Eigen::ArrayXd vSolution = ve.getPosition();
    //std::cout << "  Verlet solution: " << vSolution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in Verlet is " << exact.getError(vSolution) << std::endl;

    Eigen::ArrayXd r3Solution = r3.getPosition();
    //std::cout << "  RK3 solution: " << r3Solution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in RK3 is " << exact.getError(r3Solution) << std::endl;

    Eigen::ArrayXd r4Solution = r4.getPosition();
    //std::cout << "  RK4 solution: " << r4Solution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in RK4 is " << exact.getError(r4Solution) << std::endl;

    Eigen::ArrayXd abSolution = ab.getPosition();
    //std::cout << "  ABM4 solution: " << abSolution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in ABM4 is " << exact.getError(abSolution) << std::endl;

    Eigen::ArrayXd sSolution = rv.getPosition();
    //std::cout << "  RichardsonVerlet solution: " << sSolution.segment(0,4).transpose() << std::endl;
    //std::cout << "  Error in RichardsonVerlet is " << exact.getError(sSolution) << std::endl;

    Eigen::ArrayXd h6Solution = h6.getPosition();
    Eigen::ArrayXd h8Solution = h8.getPosition();
    Eigen::ArrayXd abbSolution = abb.getPosition();

    std::cout << maxSteps << "\t" << exact.getError(eSolution);
    std::cout << "\t" << exact.getError(r2Solution);
    std::cout << "\t" << exact.getError(a2Solution);
    std::cout << "\t" << exact.getError(vSolution);
    std::cout << "\t" << exact.getError(r3Solution);
    std::cout << "\t" << exact.getError(r4Solution);
    std::cout << "\t" << exact.getError(abSolution);
    std::cout << "\t" << exact.getError(sSolution);
    std::cout << "\t" << exact.getError(h6Solution);
    //std::cout << "\t" << exact.getError(h8Solution);
    std::cout << "\t" << exact.getError(abbSolution);
    std::cout << std::endl;
  }

  exit(0);
}
