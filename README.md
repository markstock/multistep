# multistep
Library for multistep and multistage integration of ODEs on vectors of numbers


### Description
This is a library and main routine for testing various multi-step and
multi-stage forward integrators for ODEs. It performs a simulation of 
100 bodies in three-dimensional gravitation using a wide range of time
step sizes for the following integrators: Euler, Runge-Kutta 2nd and 4th
order, Adams-Bashforth 2nd and 4th and 5th order, Standard Verlet,
a higher-order Richardson-Verlet, and a method from Hamming's "Numerical Methods for
Scientists and Engineers."

The Richardson-Verlet integrator may appear here for the first time.
It is a 4th order method using accelerations only,
manages better error than AB4 with the same history, and much better error than
RK4 when the force calculation dominates the computational effort.
It's only disadvantage is that for very small time step sizes, 
roundoff error in a subtraction prevents the total error from continuing 
to drop, as it does with AB4 and RK4. But it needs 1e-11 relative
errors for that effect to surface. It excels at RMS errors of 1e-4 to
1e-6, where it delivers similar accuracy with 1/3rd the computational
effort of RK4, and no increased error at very large time steps like
AB5.

### Compile and run
Compile and run multistep with the following commands on an RPM-based system:

    sudo yum install eigen3-devel
    g++ -o runmultistep -Ofast -std=c++11 -I/usr/include/eigen3 multistep.cpp
    ./runmultistep

### Performance
Short story: for gravitational systems, Verlet and Richardson-Verlet are the best.

### Future work
This is still a toy program, meant to test various forward integrators.
In the future, I hope to do the following:

* break it out into separate files
* use a proper build system and make a linkable library
* support multithreading and superscalar instructions for the nbody code
* clean up the code and simplify some usage
* should calculation of the highest derivative come at the end of the current step or the beginning of the next one? The former would make for cleaner code.
* include integrators with automatic time step adjustment (global first)
* develop integrator with local (element-wise) time-stepping
* calculate total effort/cost for integrators and derivative-finder
* support systems which have forcing terms on derivatives other than the highest
