
## Coding the Implicit Euler method

In every time-step we have to solve for the new value $y^{n+1}$:

$$
y^{n+1} - y^n - h f(y^{n+1}) = 0
$$

The function $f : {\mathbb R}^n \rightarrow {\mathbb R}^n$ is the right hand side of the ODE,
which has been brought into autonomous form.

We use our function algebra to build this composed function, and throw it into the Newton solver.
If we make the time-step not too large, the value $y^n$ of the old time-step is a good starting value.

To express that the independent varible is $y^{n+1}$, we create an `IdentityFunction`. The old
value is considered as a constant, i.e. a `ConstantFunction`. The right hand side function is
given by the user. Then the implicit Euler method is coded up like that:

```cpp
void SolveODE_IE(double tend, int steps,
                 VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                 std::function<void(double,VectorView<double>)> callback = nullptr)
{
  double dt = tend/steps;
  auto yold = make_shared<ConstantFunction>(y);
  auto ynew = make_shared<IdentityFunction>(y.Size());
  auto equ = ynew - yold - dt * rhs;

  double t = 0;
  for (int i = 0; i < steps; i++)
    {
      NewtonSolver (equ, y);
      yold->Set(y);
      t += dt;
      if (callback) callback(t, y);
    }
}
```

## Using the time-stepping method

A mass attached to a spring is described by the ODE

$$
m y^{\prime \prime}(t) = -k y(t)
$$

where $m$ is mass, $k$ is the stiffness of the spring, and $y(t)$ is the
displacement of the mass. The equation comes from Newton's law

> force = mass $\times$ acceleration 

We replace the second-order equation with a first order system. For mathematical simplicity we set $k = m = 1$.
Then we can define the right-hand-side as a `NonlinearFunction`. The derivative is needed for the Newton solver:

```cpp
class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};
```

Finally, We start the time-stepper with the time interval, number of steps, initial values,
right-hand-side function, and a call-back function called at the end of every time-step:
```cpp
double tend = 4*M_PI;
int steps = 100;
Vector<> y { 1, 0 };
auto rhs = make_shared<MassSpring>();

SolveODE_IE(tend, steps, y, rhs,
            [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
```
This example is provided in `demos/test_ode.cc`.


## Excercises

* Implement an explicit Euler time-stepper, and the Crank-Nicolson method.

* Compare the results of the mass-spring system for these methods, and various time-steps. Plot the solution function. What do you observe ?

* Model an electric network by an ODE. Bring it to autonomous form.
Solve the ODE numerically for various parameters with the three methods, and various time-steps.

```{image} pictures/RC.png
:width: 200px
:align: center
```


Voltage source $U_0(t) = \cos(100 \pi t)$, $R = C = 1$ or $R = 100, C = 10^{-6}$.

Ohm's law for a resistor $R$ with resistivity $R$:

$$
U = R I
$$

Equation for a capacitor $C$ with capacity $C$:

$$
I = C \frac{dU }{dt}
$$

Kirchhoff's laws:
* Currents in a node sum up to zero.
  Thus we have a constant current along the loop.
* Voltages around a loop sum up to zero. This gives:

  $$
  U_0 = U_R + U_C
  $$

Together:

$$
U_C(t) + R C \frac{dU_C}{dt}(t) = U_0(t)
$$

Use initial condition for voltage at capacitor $U_C(t_0) = 0$, for $t_0=0$.



## Stability function

An ODE $y^\prime(t) = A y(t)$ with $A \in {\mathbb R}^n$ diagonizable can be brought to $n$ scalar ODEs

$$
y_i^\prime(t) = \lambda_i y_i(t),
$$

where $\lambda_i$ are eigenvalues of $A$. If $\lambda_i$ has negative (non-positive) real part, the solution
is decaying (non-increasing). Skip index $i$.


The explicit Euler method with time-step $h$ leads to

$$
y_{i+1} = (1 + h \lambda) y_i,
$$

the implicit Euler method to

$$
y_{i+1} = \frac{1}{1-h \lambda} y_i,
$$

and the Crank-Nicolson to

$$
y_{i+1} = \frac{ 2 + h \lambda }  { 2 - h \lambda } y_i
$$

The stability function $g(\cdot)$ of a method is defined such that

$$
y_{i+1} = g(h \lambda) y_i
$$

These are for the explicit Euler, the implicit Euler, and the Crank-Nicolson:

$$
g_{EE}(z) = 1+z \qquad g_{IE}(z) = \frac{1}{1-z}  \qquad g_{CN}(z) = \frac{2+z}{2-z}
$$

The domain of stability is

$$
S = \{ z : | g(z) | \leq 1 \}
$$

What is

$$
S_{EE} = \{ z : |z + 1| \leq 1 \}
$$

$$
S_{IE} = \{ z : | z - 1 | \geq 1 \}
$$

$$
S_{CN} = \{ z : {\text Re} (z) \leq 0 \}
$$

If Re$(\lambda) \leq 0$, then $h \lambda$ is always in the domain of stability of the implicit Euler,
and of the Crank-Nicolson. This property of an method is called $A$-stability. The explicit Euler leads to
(quickly) increasing numerical solutions if $h$ is not small enough.

If $\lim_{z \rightarrow -\infty} g(z) = 0$, quickly decaying solutions lead to quickly decaying
numerical solutions. This is the case for the implicit Euler, but not for the Crank-Nicolson. This property is called $L$-stability. One observes (slowly decreasing) oscillations with the CR - method when $-h \lambda$ is large.

