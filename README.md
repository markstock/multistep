# multistep
Library for multistep and multistage integration of ODEs on vectors of numbers


### Description
This is a library and main routine for testing various multi-step and
multi-stage forward integrators for ODEs. It performs a simulation of 
100 bodies in three-dimensional gravitation using a wide range of time
step sizes for the following integrators: Euler, Runge-Kutta 2nd and 4th
order, Adams-Bashforth 2nd and 4th and 5th order, Standard Verlet,
higher-order Verlet, and a method from Hamming's "Numerical Methods for
Scientists and Engineers."

### Compile and run
Compile and run multistep with the following commands on an RPM-based system:

    sudo yum install eigen3-devel
    g++ -o runmultistep -Ofast -std=c++11 -I/usr/include/eigen3 multistep.cpp
    ./runmultistep

### Performance
Coming soon
