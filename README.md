Physics LIB<a name="TOP"></a><br>[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Computational-science-based-on-HPC/Physics_LIB/blob/master/LICENSE)
=============================

This Repository Contains **Physics LIB**, a physics library that simulates various physics calculations as _parallelelized_ and _seriallized_ implementations to obtain optimal performance and results.
- - - - 
# Features #
### 1- Damped spring motion ###
  - _Damped Spring Motion_ is a __simple harmonic motion__ that simulates the motion of a mass attached to a spring and oscillates in a y-axis direction with given parameters that the user specifies like mass and stiffness of spring.
  - This simulation is implemented in both parallel and serial implementations.
### 2- Elastic Pendulum ###
  - _Elastic pendulum/pendulum spring / swinging spring system_ is a system where it simulates the motion of 2D-spring in both directions y-axis and x-axis = given parameters that the user specifies like mass and stiffness of spring.
  - Due to the chaotic motion of this system and its dependencies its algorithm is implemented in serial algorithm only.
### 3- 1D Heat Equation ###
  - Heat Equation is an equation where it simulates the change of heat in the body isolated from the outer world and isn't affected by any other heat or cooling source.
  - In this simulation, we simulate heat propagation in a 1D body as _wire_ and the change in its temperature over time using Fourier transform.
  - This simulation is implemented as a serial and parallel program.
### 4- 2D Heat Equation ###
  - Heat Equation is an equation where it simulates the change of heat in the body isolated from the outer world and isn't affected by any other heat or cooling source.
  - In this simulation, we simulate heat propagation in a 2D body as _cpu_ as square or rectangle and the change in its temperature over time using Fourier transform.
  - This simulation is implemented as a serial and parallel program.
- - - - 
# Technologies and Tools Used #
  - Programming Language: C
  - Parallel Computing: MPI, OpenMP
  - Build: Make
- - - - 
# Dynamic Library #
## Installation ## 
 ```  
make
sudo make install
make clean
```
To include in the environment path:
```
sudo ldconfig
echo $LD_LIBRARY_PATH
```
If couldn't find the library location listed then create a new environment variable:
```
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATHD
```
Now it has been added to environment variables to be sure:
```
echo $LD_LIBRARY_PATH
```
Make sure it has been listed.
## Compilation and Run ##
 ``` 
mpicc -o myprogram main.c -lphysics -fopenmp
mpiexec -n 4 ./myprogram
```
- - - - 
# Static Library #
## Installation ## 
 ``` 
cd ./make 
make
make clean
```
## Compilation and Run ##
 ``` 
 mpicc ../examples/main.c -o main -L../make -l:libphysics.a -fopenmp -lm -lc
 mpiexec -n 4 ./main
```
- - - - 
## Heat Equation 1D ##
The following graph represents the execution time of the _Heat Equation 1D_ algorithm via different implementations (MPI, OpenMP, Serial) and the change of it when changing the **Precision** _(number of iterations)_.


![Heat Equation 1D](https://raw.githubusercontent.com/Computational-science-based-on-HPC/Physics_LIB/master/graphs/heat1d%20Graph.png)
## Heat Equation 2D ##
The following graph represents the execution time of the _Heat Equation 2D_ algorithm via different implementations (MPI, OpenMP, Serial) and the change of it when changing the **Precision** _(number of iterations)_.


![Heat Equation 2D graph](https://raw.githubusercontent.com/Computational-science-based-on-HPC/Physics_LIB/master/graphs/Heat%202d%20Graph.png)
## Damped Oscillation ##
The following graph represents the execution time of the _Damped Oscillation Simulation_ algorithm via different implementations (MPI, OpenMP, Serial) and the change of it when changing the **Time Limit** _(number of iterations)_, and discusses the execution time of the simulation when false sharing occurred.


![Damped Oscillation Graph](https://raw.githubusercontent.com/Computational-science-based-on-HPC/Physics_LIB/master/graphs/DampedOSC%20Graph.png)
- - - -
# Contributing
We welcome contributions from the community to improve our project. Please read our [Contributing Guidelines](https://github.com/Computational-science-based-on-HPC/Physics_LIB/blob/master/CONTRIBUTING.md) to understand how you can contribute and make the process smoother for everyone involved.

# Code of Conduct

We have a [Code of Conduct](https://github.com/Computational-science-based-on-HPC/Physics_LIB/blob/master/CODE_OF_CONDUCT.md) in place to promote a positive and inclusive community. We expect all contributors to adhere to the code when participating in discussions or contributing to the project.

# License #
This project is licensed under the [MIT License](https://github.com/Computational-science-based-on-HPC/Physics_LIB/blob/master/LICENSE).

## Links ## 
- [Documentation](https://Computational-science-based-on-HPC.github.io/index.html)
- [GitHub](https://github.com/Computational-science-based-on-HPC/Physics_LIB)

