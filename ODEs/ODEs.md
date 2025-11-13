# Solving ordinary differential equations

We are developing a library for solving nonlinear ordinary differential equations. First we build an abstraction of nonlinear functions, and implement solvers for nonlinear equations. On this basis we can simply define time-stepping algorithms. We apply it to multi-body systems. OpenGL visualization is provided based on the python interface to threejs.

- Solving non-linear equations
- Explicit and implicit time-stepping methods
- Connecting masses with springs
- Introducing constraints



Fork that github repo [ASC-ODE](https://github.com/TUWien-ASC/ASC-ODE.git), and clone it via

    git clone https://github.com/my_fork_on_github/ASC-ODE.git

It contains the basic linear algebra package nanoblas as an git submodule. Populate it using

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

