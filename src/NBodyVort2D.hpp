/*
 * NBodyVort2D.hpp - a simple 2D vortex method
 *
 * Copyright 2016,22 Mark J. Stock, markjstock@gmail.com
 */

#include "DynamicalSystem.hpp"

#include <Eigen/Dense>

#include <cstdint>


/*
 * A 2D vortex system with constant strengths, radii
 */
class NBodyVort2D : public VelocitySystem<Eigen::ArrayXd> {
public:
  NBodyVort2D(const int32_t _num) :
    VelocitySystem<Eigen::ArrayXd>(2*_num),
    num(_num)
  {
    // num is how many bodies

    // we store the unchanging properties here
    circ = Eigen::ArrayXd::Random(num) / (2.0*num);
    radiusSquared = 0.2 + 0.1 * Eigen::ArrayXd::Random(num);
    radiusSquared = radiusSquared.square();

    // store initial conditions in case we want to reuse this
    ic.x[0] = 10.0 * Eigen::ArrayXd::Random(numVars);
  };

  // perform n-body acceleration calculation; uses position and mass and radius squared
  Eigen::ArrayXd getHighestDeriv(const Eigen::ArrayXd pos) {

    // create the accumulator vector
    Eigen::ArrayXd newVels = Eigen::ArrayXd::Zero(numVars);

    // evaluate the new vels
    if (true) {
      for (int32_t i=0; i<num; ++i) {
        // new velocities on particle i
        Eigen::Vector2d thisVel(0.0, 0.0);
        for (int32_t j=0; j<num; ++j) {
          // 20 flops
          // the influence of particle j
          Eigen::Vector2d dx = pos.segment(2*j,2) - pos.segment(2*i,2);
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
  int32_t num;
  // these values are constant in time and unique to this system
  Eigen::ArrayXd circ;
  Eigen::ArrayXd radiusSquared;
};

