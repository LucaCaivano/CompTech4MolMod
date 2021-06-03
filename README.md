# Computational Techniques for Molecular Modeling
Projects of the numerical part of the course *Computational Techniques for Molecular Modeling* @ PoliMi - First Edition, Spring 2021

## Team

* Luca Caivano (Mathematical Engineering)
* Andrea Della Libera (Chemical Engineering)
* Manfred Nesti (Mathematical Engineering)
* Valeria Pajola (Chemical Engineering)
* Alessandro Pegurri (Chemical Engineering)
* Bruno Ursino (Mathematical Engineering)
* Chiara Vitale (Chemical Engineering)

## Assignments

### Assignment 1
* Implement a simulator of the conservative double pendulum using the Velocity-Verlet method as symplecting integrator and study the accuracy of trajectories and total energy with respect to time step.
* Reformulate the model so that the contraints are accounted for, simulating in cartesian coordinates using the SHAKE method.
* Apply the Velocity-Verlet method to the Frozen Argon Crystal example as in *Hairer, Lubich, Wanner - Geometric Numerical Integration - Section I.4*.

### Assignment 2
Reproduce the example in *Gabriel, Knapek, Zumbusch - Numerical Simulation in Molecular Dynamics - Section 3.6.1* using periodic boundary conditions.

### Assignment 3
Reproduce the example of Rayleigh-Taylor instability with Coulomb potential in *Gabriel, Knapek, Zumbusch - Numerical Simulation in Molecular Dynamics - Section 7.4.1*.

### Assignment 4
For the C2H6 molecule

* Compute the equilibrium using GAUSSIAN
* Compute forces by invoking GAUSSIAN and parse output via Octave / Bash script
* Implement a descent method in Octave to perform the optimization using for example Gradient, BFG and Nonlinear Conjugate Gradien as descent directions
* Apply box constraints via orthogonal projection with feasible Projected Gradient or Projected quasi-Newton
* Apply Reduced Rank Extrapolation

#### Assignment 4.5 _(only for 5CFU course version)_
Implement also the Quadratic Line Search
([video reference](https://youtu.be/MKmIvtq83LY))
