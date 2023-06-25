libphysics<a name="TOP"></a>
===================

- - - - 
## Description ##
**libphysics** is a physics library that simulate various physics calculation as _parallelelized_ and _seriallized_ implementations to obtain optimal performance and results.
- - - - 
## Simulations ##
### 1- Damped spring motion ###
  - _Damped Spring Motion_ is a __simple harmonic motion__ that simulates the motion of mass attached to spring and oscillates in y-axis direction with given parameters that the user specify like mass and stiffness of spring.
  - This simulation is implemented in both parallel and serial implementations.
### 2- Elastic Pendulum ###
  - _Elastic pendulum / pendulum spring / swinging spring system_ is a system where it simulates the motion of 2D-spring in both direction y-axis and x-axis = given parameters that the user specify like mass and stiffness of spring .
  - Due to the chaotic motion of this system and dependencies its algorithm is implemented in serial algorithm only.
### 3- 1D Heat Equation ###
  - Heat Equation is equation where it simulate the change of heat in body isolated from outer world and doesn't affected by any other heat or cooling source.
  - In this simulation we simulate heat propagation in 1D body as _wire_ and the change in its tempreture over time using fourier transform.
  - This simulation is implemented as serial and parallel program.
### 4- 2D Heat Equation ###
  - Heat Equation is equation where it simulate the change of heat in body isolated from outer world and doesn't affected by any other heat or cooling source.
  - In this simulation we simulate heat propagation in 2D body as _cpu_ as square or rectangle and the change in its tempreture over time using fourier transform.
  - This simulation is implemented as serial and parallel program.
  - - - - - 
