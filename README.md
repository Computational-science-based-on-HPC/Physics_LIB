Physics LIB<a name="TOP"></a>
===================

- - - - 
## Description ##
**Physics LIB** is a physics library that simulates various physics calculations as _parallelelized_ and _seriallized_ implementations to obtain optimal performance and results.
- - - - 
## Simulations ##
### 1- Damped spring motion ###
  - _Damped Spring Motion_ is a __simple harmonic motion__ that simulates the motion of a mass attached to a spring and oscillates in a y-axis direction with given parameters that the user specifies like mass and stiffness of spring.
  - This simulation is implemented in both parallel and serial implementations.
### 2- Elastic Pendulum ###
  - _Elastic pendulum/pendulum spring / swinging spring system_ is a system where it simulates the motion of 2D-spring in both directions y-axis and x-axis = given parameters that the user specifies like mass and stiffness of spring.
  - Due to the chaotic motion of this system and its dependencies its algorithm is implemented in serial algorithm only.
### 3- 1D Heat Equation ###
  - Heat Equation is an equation where it simulates the change of heat in body isolated from outer world and isn't affected by any other heat or cooling source.
  - In this simulation, we simulate heat propagation in a 1D body as _wire_ and the change in its temperature over time using Fourier transform.
  - This simulation is implemented as a serial and parallel program.
### 4- 2D Heat Equation ###
  - Heat Equation is an equation where it simulates the change of heat in body isolated from outer world and isn't affected by any other heat or cooling source.
  - In this simulation, we simulate heat propagation in a 2D body as _cpu_ as square or rectangle and the change in its tempreture over time using Fourier transform.
  - This simulation is implemented as a serial and parallel program.
- - - - 
# Dynamic Library #
## Installation ## 
 ```  
make
make install
make clean
```
To include in environment path:
```
sudo ldconfig
echo $LD_LIBRARY_PATH
```
If couldn't find the library location listed then create new environment variable:
```
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATHD
```
Now it been added to environment variables to be sure:
```
echo $LD_LIBRARY_PATH
```
Make sure it been listed.
## Compilation and Run ##
 ``` 
mpicc -o myprogram main.c -L/usr/local/lib/ -lphysics -fopenmp
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
## Links ## 
- [Documentation](https://Computational-science-based-on-HPC.github.io/index.html)
- [GitHub](https://github.com/Computational-science-based-on-HPC/Physics_LIB)
