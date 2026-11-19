# Implementing time-stepping methods

We implement the methods from the previous section in the C++ library ASC-ODE.
The library is built from three layers:

1. **Functions**: the class `NonlinearFunction` represents a function $f : {\mathbb R}^{m} \rightarrow {\mathbb R}^{n}$
   which can evaluate its value $f(x)$ and its Jacobi-matrix $f^\prime(x)$.
   The right-hand side of an ODE is such a function.
   Simple building blocks (identity, constants, sums, scaling, composition) are combined to
   more complex functions. We call this the function algebra.

2. **Solvers**: Newton's method solves $F(x) = 0$ for any `NonlinearFunction` $F$.
   It needs nothing else than the two evaluation functions.

3. **Time-steppers**: the class `TimeStepper` has the virtual function `doStep(tau, y)`,
   which advances the state `y` by one time-step `tau`.
   An explicit method just evaluates the right-hand side. An implicit method sets up
   the equation for the new state as an expression of the function algebra, and
   calls Newton's method.

The essential point is how these layers work together: For the implicit Euler method, the new state is the solution of

$$
F(y) := y - y_{\text{old}} - \tau f(y) = 0.
$$

Newton's method needs the Jacobi-matrix $F^\prime(y) = I - \tau f^\prime(y)$. We never implement this formula.
The expression `ynew - yold - tau * rhs` builds $F$ from building blocks, and every building block
knows how to differentiate itself by the sum-, product-, and chain-rules. The same mechanism
gives us Crank-Nicolson, and later implicit Runge-Kutta methods, without additional derivative code.
The user only provides $f$ and $f^\prime$, and in {doc}`implementation_ad` we will see how
to obtain $f^\prime$ automatically.

The code is organized in these files:

| file | content |
|------|---------|
| `src/nonlinfunc.hpp` | `NonlinearFunction` and the building blocks of the function algebra |
| `src/Newton.hpp` | `NewtonSolver` |
| `src/timestepper.hpp` | `TimeStepper`, `ExplicitEuler`, `ImplicitEuler` |
| `src/autodiff.hpp` | the `AutoDiff` class for automatic differentiation |
| `demos/test_ode.cpp` | the mass-spring example |
| `demos/plotmassspring.py` | plotting the results |

The following sections discuss these parts:

```{tableofcontents}
```
