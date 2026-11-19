# Mechanical Systems

We consider a system of point masses $m_i$ at time-dependent positions $x_i(t)$.

There is potential energy $U(x)$ stored in that system. Examples are gravitational energy

$$
\sum m_i g x_{i,z},
$$

or energy stored in springs $k_{ij}$ connecting masses $m_i$ and $m_j$:

$$
\sum_{k_{ij}} \frac{1}{2} k_{ij} ( |x_i - x_j | - l_{ij} )^2,
$$

where $l_{ij}$ is the length of the spring without load, and $|x_i - x_j|$ is the Euclidean distance.


Taking the negative derivative

$$
F_i := -\frac{\partial U}{\partial x_i}(x)
$$

gives the force acting on the mass point $i$.

In case of gravity it is $-m_i g e_z$, in case of springs it is

$$
-\sum_j k_{ij} ( |x_i - x_j | - l_{ij} ) \frac{x_i - x_j}{|x_i - x_j|}
$$

i.e. spring constant multiplied by elongation, in the direction of the vector $x_j-x_i$.
The sum runs over all mass points $j$ connected to mass point $i$ with a spring.


By Newton's principle, the force leads to acceleration:

$$
m_i \ddot x_i = -\frac{\partial U}{\partial x_i}
$$

Here, $\dot {x}$ denotes the derivative w.r.t. time $t$.

## Energy conservation

The total energy is the sum of kinetic and potential energy:

$$
E(t) = \sum_i \frac{1}{2} m_i |\dot x_i(t)|^2 + U(x(t)).
$$

Along a solution it is constant:

$$
\frac{d}{dt} E = \sum_i m_i \, \dot x_i \cdot \ddot x_i + \sum_i \frac{\partial U}{\partial x_i} \cdot \dot x_i
= \sum_i \Big( m_i \ddot x_i + \frac{\partial U}{\partial x_i} \Big) \cdot \dot x_i = 0.
$$

For simulations over long times, a numerical method should not destroy this property.
We have seen in section {ref}`sec:stability` for the oscillator $y^{\prime \prime} = -\omega^2 y$ that the explicit Euler method increases the
energy, the implicit Euler method damps it, and the Crank-Nicolson method preserves it exactly.
For non-linear problems, so called symplectic methods (for example the implicit mid-point rule, Gauss-Legendre methods, and
the Störmer-Verlet method below) do not preserve the energy exactly. But the energy error stays bounded by $O(\tau^p)$ over very long times,
there is no drift.

## The Newmark method

We start with reducing the second order equation to a first order system. Therefore, we introduce velocity $v_i = \dot x_i$:

$$
\begin{array}{rcl}
\dot x_i & = & v_i \\
m_i \dot v_i & = & F_i(x)
\end{array}
$$

Applying the Crank-Nicolson time-stepping scheme we obtain

$$
\begin{array}{rcl}
x_{n+1} - x_n & = & \frac{\tau}{2} (v_n + v_{n+1} ) \\
M (v_{n+1} - v_n) & = & \frac{\tau}{2} (F(x_n) + F(x_{n+1}))
\end{array}
$$

By renaming $a_n := M^{-1} F(x_n)$, we can substitute 

$$
\begin{array}{rcl}
v_{n+1} & = & v_n + \frac{\tau}{2} (a_n + a_{n+1}) \\
x_{n+1} & = & x_n + \frac{\tau}{2} (v_n + v_{n+1}) = x_n + \tau v_n + \frac{\tau^2}{4} (a_n + a_{n+1} )
\end{array}
$$

and obtain the equation for the new $a_{n+1}$:

$$
M a_{n+1} = F \big( x_n + \tau v_n + \tfrac{\tau^2}{4} (a_n + a_{n+1} ) \big)
$$

The dimension of the nonlinear system is the same as the number of state variables (and not twice the size,
as it looks after introducing the first order system).

The function `SolveODE_Newmark` in `mechsystem/Newmark.hpp` (branch `Newmark`) solves $M \ddot x = F(x)$. Both $M$ and $F$
are given as `NonlinearFunction`s. The new acceleration is the independent variable, and the function algebra
expresses the new velocity and position as functions of it (with $\gamma = \frac{1}{2}$, $\beta = \frac{1}{4}$):

```cpp
auto anew = std::make_shared<IdentityFunction>(a.size());
auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);

auto equ = Compose(mass, anew) - Compose(rhs, xnew);
```

In every time-step, Newton's method solves `equ` for the new acceleration. The same function objects `xnew` and `vnew` then
evaluate the new position and velocity:

```cpp
NewtonSolver (equ, a);
xnew -> evaluate (a, x);
vnew -> evaluate (a, v);

xold->set(x);
vold->set(v);
aold->set(a);
```

## The Störmer-Verlet method

An explicit alternative is the Störmer-Verlet method:

$$
\begin{array}{rcl}
v_{n+1/2} & = & v_n + \frac{\tau}{2} a(x_n) \\
x_{n+1} & = & x_n + \tau \, v_{n+1/2} \\
v_{n+1} & = & v_{n+1/2} + \frac{\tau}{2} a(x_{n+1})
\end{array}
$$

with $a(x) = M^{-1} F(x)$. In the notation of the generalized-$\alpha$ method below, it is the Newmark family with $\beta = 0$, $\gamma = \frac{1}{2}$,
i.e. $x_{n+1} = x_n + \tau v_n + \frac{\tau^2}{2} a_n$, while the method above has $\beta = \frac{1}{4}$, $\gamma = \frac{1}{2}$.
It needs one evaluation of the force per time-step, since $a(x_{n+1})$ is reused in the next step.
The method is of second order, symplectic, and time-reversible. Eliminating the velocities gives the two-step form

$$
x_{n+1} - 2 x_n + x_{n-1} = \tau^2 a(x_n),
$$

which we use for the wave equation in {doc}`../PDEs/timedependent/verlet`.
As an explicit method it is only conditionally stable: for the oscillator $\ddot x = -\omega^2 x$ one needs $\tau \omega < 2$.
Stiff springs therefore require small time-steps, and implicit methods are preferred.

## Generalized-$\alpha$ method

Although the Newmark method preserves energy exactly for linear ODEs, it leads to
instabilities for non-linear equations. For this reason, the following generalization is
usually used in mechanical simulation codes:

By adding a bunch of parameters, one ends up with the *Generalized $\alpha$ method*:

$$
M a_{n+1-\alpha_m} = F(x_{n+1-\alpha_f})
$$

with

$$
\begin{array}{rcl}
x_{n+1} & = & x_n + \tau v_n  + \tau^2 ( ( \frac{1}{2} - \beta) a_n + \beta a_{n+1} ) \\
v_{n+1} & = & v_n + \tau ( (1-\gamma) a_n + \gamma a_{n+1} ) \\
x_{n+1-\alpha_f} & = & (1-\alpha_f) x_{n+1} + \alpha_f x_n \\
a_{n+1-\alpha_m} & = & (1-\alpha_m) a_{n+1} + \alpha_m a_n 
\end{array}
$$

Depending on one free parameter $\rho^\infty$, parameters tuned for optimal accuracy
and stability are obtained by

$$
\begin{array}{rcl}
\alpha_m & = & \frac{2 \rho^\infty - 1}{\rho^\infty + 1} \\
\alpha_f & = & \frac{\rho^\infty}{\rho^\infty + 1} \\
\beta & = & \frac{1}{4} (1 - \alpha_m + \alpha_f)^2 \\
\gamma & = & \frac{1}{2} - \alpha_m + \alpha_f
\end{array}
$$

The parameter $\rho^\infty \in [0,1]$ specifies the damping for high frequency functions: $\rho^\infty = 1$
gives no numerical damping (like the Newmark method), smaller values damp high frequencies stronger.
The implementation `SolveODE_Alpha` in ASC-ODE uses the variant $(1-\alpha_f) F(x_{n+1}) + \alpha_f F(x_n)$
instead of $F(x_{n+1-\alpha_f})$, which coincides for linear $F$.
The reference is J. Chung and G.M. Hulbert *A Time Integration Algorithm for Structural Dynamics With Improved Numerical Dissipation: The Generalized-α Method*, J. Appl. Mech., 1993, [see](https://asmedigitalcollection.asme.org/appliedmechanics/article/60/2/371/423023/A-Time-Integration-Algorithm-for-Structural?casa_token=PkMuU1hl_4AAAAAA:btyl_IyLpSdrI9q5LjNYSDLzOECUIoT1b3IRqcF86A773ZNWRnyuN94Tizg4kQUaVw2gZftWhg)

An easy to read derivation is [here](https://miaodi.github.io/finite%20element%20method/newmark-generalized/).




## Systems with constraints

Think of a pendulum, where a mass $m$ is attached to a (mass-less) stick of length $l$.
One can describe the system by the angle $\alpha$ between the stick and the vertical axis.
By some high-school mathematics, one splits the total gravitational force into the direction of
the stick, and the orthogonal complement, and ends up with the equation for the angle $\alpha$:

$$
\ddot \alpha = -\frac{g}{l} \sin \alpha
$$

```{figure} pictures/pendulum.png
:width: 55%
:align: center

Pendulum with anchor point $x_0$, stick of length $l$, and angle $\alpha$. The constraint $c(x) = 0$ keeps the
mass on a circle. Gravity $-m g\, e_z$ and the constraint force $\lambda \nabla_x c$ act on the mass.
```

However, such a direct modeling with angles (so called generalized coordinates) becomes
difficult for more complex systems. A simpler possibility is to pose the length constraint

$$
0 = c(x) := | x - x_0 |^2 - l^2,
$$

where $x_0$ is the anchor point of the pendulum, and define the Lagrange function

$$
L(x, \lambda) = -U(x) + \left< \lambda , c(x) \right>,
$$

with $U(x) = m g x_z$, and pose the equation

$$
\begin{array}{rcl}
m_i \ddot x_i & = & \frac{\partial}{\partial x_i} L(x, \lambda) \\
 0 & = & \nabla_\lambda L(x, \lambda)
\end{array}
$$

The physical meaning of the Lagrange parameter $\lambda$ is the longitudinal force in the pendulum (precisely: $\lambda \nabla_x c$ is the force vector).
This is an ordinary differential equation, with algebraic constraints,
a so called differential-algebraic equation (DAE).

One can directly use the generalized-$\alpha$ method for this DAE. The right-hand side is the gradient of the Lagrange function, the mass matrix on the left-hand side is extended by $0$ for the Lagrange parameters:

$$
\left( \begin{array}{cc}
 M & 0 \\
 0 & 0
 \end{array} \right)
\left( \begin{array}{c}
   \ddot x \\
   0
 \end{array} \right) =
 \nabla L (x, \lambda)
$$

The Lagrange parameter is treated like an additional position variable. Its "acceleration" is an auxiliary unknown,
only the resulting $\lambda_{n+1}$ has a meaning. With the variant used in ASC-ODE, the rows of the Lagrange parameters enforce

$$
(1-\alpha_f) \, c(x_{n+1}) + \alpha_f \, c(x_n) = 0,
$$

thus the constraint is satisfied exactly in every time-step, provided it is satisfied initially.

Some practical consequences of the missing equation for $\dot \lambda$:

* Explicit methods cannot be applied directly.

* The initial values must be consistent: $c(x(0)) = 0$, and also the time derivative $\nabla c(x(0)) \cdot \dot x(0) = 0$,
  i.e. the initial velocity is tangential to the constraint.

* The Lagrange parameter is sensitive to perturbations. Without numerical damping ($\rho^\infty = 1$), $\lambda$ oscillates
  from time-step to time-step, while the positions are fine. Numerical damping with $\rho^\infty < 1$ removes these oscillations.


## Implementing mass-spring systems

Merge branch `Newmark` from ASC-ODE into your main branch, see [here](https://github.com/TUWien-ASC/ASC-ODE/tree/Newmark).

For small systems it is reasonable to implement the right-hand-side function directly.
For larger mechanical systems we describe the system by data structures, which are
implemented in `mechsystem/mass_spring.hpp`. They are templated on the spatial dimension `D`:

```cpp
template <int D> class Mass { public: double mass; Vec<D> pos; Vec<D> vel = 0.0; Vec<D> acc = 0.0; };
template <int D> class Fix  { public: Vec<D> pos; };

class Connector { public: enum CONTYPE { FIX=1, MASS=2 }; CONTYPE type; size_t nr; };
class Spring { public: double length; double stiffness; std::array<Connector,2> connectors; };
```

A `Connector` refers to either a fixed point or a mass, by its type and index. The class `MassSpringSystem<D>`
has containers for fixes, masses and springs, and the gravity vector. The functions `addFix` and `addMass` return
connectors, which are then used to add springs.

The state vector of the ODE consists of the positions of all masses, $D$ values per mass.
The class `MSS_Function<D>` evaluates the right-hand side for the `MassSpringSystem`. It interprets the state vector as a matrix
with one row per mass (via `asMatrix`), starts with the gravitational forces, and adds the forces of all springs:

```cpp
for (auto spring : mss.springs())
  {
    auto [c1,c2] = spring.connectors;
    Vec<D> p1, p2;
    if (c1.type == Connector::FIX)
      p1 = mss.fixes()[c1.nr].pos;
    else
      p1 = xmat.row(c1.nr);
    if (c2.type == Connector::FIX)
      p2 = mss.fixes()[c2.nr].pos;
    else
      p2 = xmat.row(c2.nr);

    double force = spring.stiffness * (norm(p1-p2)-spring.length);
    Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
    if (c1.type == Connector::MASS)
      fmat.row(c1.nr) += force*dir12;
    if (c2.type == Connector::MASS)
      fmat.row(c2.nr) -= force*dir12;
  }
```

Finally every row is divided by its mass. Thus `MSS_Function` returns the acceleration $M^{-1} F(x)$, and the mass
function passed to the time-stepping method is the `IdentityFunction`. The derivative is computed by numerical differentiation.

Python bindings in `bind_mass_spring.cpp` allow to set up the system from Python. The function `simulate(tend, steps)`
calls the generalized-$\alpha$ method with $\rho^\infty = 0.8$, and stores the new state in the masses:

```python
from mass_spring import *

mss = MassSpringSystem3d()
mss.gravity = (0,0,-9.81)

mA = mss.add (Mass(1, (1,0,0)))
mB = mss.add (Mass(2, (2,0,0)))
f1 = mss.add (Fix( (0,0,0)) )
mss.add (Spring(1, 10, (f1, mA)))
mss.add (Spring(1, 20, (mA, mB)))

mss.simulate (0.1, 10)
print ("state = ", mss.getState())
```

The directory `mechsystem` contains
- file `mass_spring.cpp` for testing a system from C++
- file `test_mass_spring.py` for testing the system from Python
- Jupyter-notebook file `mass_spring.ipynb` for testing the system including 3D visualization


## Extending the code

The following extensions are your own implementation, here are suggestions for the interfaces.

**Distance constraints.** A constraint fixes the distance between two masses, or between a mass and a fixed point:

```cpp
class DistanceConstraint { public: double length; std::array<Connector,2> connectors; };
```

The `MassSpringSystem` gets a container of constraints and a function `addConstraint`. The state vector is extended by
one Lagrange parameter per constraint, its dimension becomes $D \cdot \#\text{masses} + \#\text{constraints}$.
For the mass rows, `MSS_Function` adds the constraint forces $\lambda \nabla_x c$ (divided by the mass, as all forces),
for the constraint rows it returns the constraint $c(x)$.
The mass function is not the identity anymore, but must vanish for the Lagrange parameters.
This is the existing `Projector`:

```cpp
auto mass = std::make_shared<Projector> (dim, 0, D*mss.masses().size());
```

**Exact derivative.** Every spring contributes $D \times D$ blocks to four positions of the Jacobi matrix of `MSS_Function`.
Test your implementation by comparing with the numerical derivative.

**Störmer-Verlet.** A function with the same arguments as `SolveODE_Newmark`, but without the mass function:

```cpp
void SolveODE_Verlet (double tend, int steps,
                      VectorView<double> x, VectorView<double> dx,
                      std::shared_ptr<NonlinearFunction> rhs,
                      std::function<void(double,VectorView<double>)> callback = nullptr);
```

Here, `rhs` returns the acceleration $M^{-1} F(x)$, as `MSS_Function` does. As an explicit method, it needs neither Newton's method nor the derivative of `rhs`.
Alternatively, design a class hierarchy of time-steppers for second order systems, with a virtual function `doStep(tau, x, v)`.
Where do you store the acceleration $a(x_{n+1})$, which is reused in the next time-step?


## Exercises

<!--
### Test Newmark and Generalized-$\alpha$ solvers

- compile and run the cods test_newmark and test_alpha. Plot the solution functions

- test the equation for the pendulum using the angle $\alpha$

- implement a double pendulum as DAE
-->

### Mass-spring system

- build and run the examples in the `mechsystem` directory
- add distance constraints to the MassSpring system. Enforce them using Lagrange multipliers. <!--, similar as in `test_alpha.cc`.-->
- the class `MSS_Function` implements the derivative using numerical differentiation. Implement the exact derivative.
- Experiment with different mechanical structures:
  * Extend the double-pendulum to a chain
  * Build more complex structures using beams (for example a crane), and simulate vibration
  * Replace the springs of the double-pendulum by distance-constraints
  * Build a spinning top (German: Kreisel). We know angular momentum keeps it spinning.
    Connect three masses by distance constraints modeling the spinning top.
    Start with velocities consistent to a rotating motion

### Energy and constraints

- Simulate the pendulum $\ddot \alpha = -\frac{g}{l} \sin \alpha$ over a long time (e.g. 1000 periods) with the explicit Euler,
  implicit Euler, Crank-Nicolson, and Störmer-Verlet methods. Plot the energy over time.
- Implement the Störmer-Verlet method for the `MassSpringSystem`. Find the largest stable time-step for a chain of masses
  and springs, and compare with the eigenvalues of the linearized system.
- Simulate the pendulum as DAE with the generalized-$\alpha$ method for $\rho^\infty = 1$ and $\rho^\infty = 0.8$.
  Plot the constraint $c(x)$ and the Lagrange parameter $\lambda$ over time.









