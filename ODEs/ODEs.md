# Solving ordinary differential equations

We are developing a library for solving nonlinear ordinary differential equations, and learn the numerical analysis behind the methods.
First we build an abstraction of nonlinear functions, and implement solvers for nonlinear equations. On this basis we can simply define time-stepping algorithms.
We apply it to mechanical systems of masses connected by springs and constraints. 3D visualization in Jupyter notebooks is provided by pythreejs.

The chapter covers:

- {doc}`theory`: existence, uniqueness, and stability of solutions
- {doc}`firstmethods`: the explicit and implicit Euler and the Crank-Nicolson methods, their convergence and stability, stiff ODEs
- {doc}`implementation`: function algebra, Newton's method, explicit and implicit time-steppers in C++
- {doc}`implementation_ad`: numerical and automatic differentiation
- {doc}`RungeKutta`: explicit and implicit Runge-Kutta methods of arbitrary order
- {doc}`adaptive`: embedded methods, and automatic step-size control for non-stiff and stiff problems
- {doc}`mechanical`: energy conservation, the Newmark, Störmer-Verlet, and generalized-$\alpha$ methods, systems with constraints


## Getting the code

Fork that github repo [ASC-ODE](https://github.com/TUWien-ASC/ASC-ODE.git), and clone it via

    git clone https://github.com/my_fork_on_github/ASC-ODE.git

It contains the basic linear algebra package nanoblas as a git submodule. Populate it using

    cd ASC-ODE
    git submodule update --init

You are welcome to replace nanoblas by your preferred linear algebra package.

ASC-ODE has the following dependencies which must be installed:

  * pybind11
  * lapack

Then you can use cmake to build ASC-ODE:

    mkdir build
    cd build
    cmake ..
    make

The Runge-Kutta methods and the mechanical systems are provided in the branches `RungeKutta` and `Newmark`.
The corresponding sections explain how to merge them into your fork.
