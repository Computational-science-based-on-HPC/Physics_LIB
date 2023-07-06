Physics LIB<a name="TOP"></a>
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
# Contribution # 
Contributions to the PhysicsLIB project are welcome! If you find any issues or have suggestions for improvements, feel free to open an issue or submit a pull request.
Please ensure that your contributions adhere to the Contributor Covenant and that you follow the project's style guide.
- - - - 
# License #
This project is licensed under the [MIT License](https://github.com/Computational-science-based-on-HPC/Physics_LIB/blob/master/LICENSE).
## Links ## 
- [Documentation](https://Computational-science-based-on-HPC.github.io/index.html)
- [GitHub](https://github.com/Computational-science-based-on-HPC/Physics_LIB)

