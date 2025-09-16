# Implementation

We implement a function algebra, which allows us to write expressions like
```cpp
  func = Id + 3 * Compose (f, g);
```
where `f` and `g` are functions, and `Id` is the identic function. Then
the composed function `func` shall be able to compute the function value and the
derivative at a given point:
```cpp
  Matrix jacobi(func->DimF(), func->DimX());
  func->EvaluateDeriv(x, jacobi);
```

The base class for such functions is
```cpp
class NonlinearFunction {
public:
  virtual ~NonlinearFunction() = default;
  virtual size_t DimX() const = 0;
  virtual size_t DimF() const = 0;
  virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
};
```
`DimX` and `DimF` provide the vector space dimensions of the domain, and the image.
The `Evaluate` and `EvaluateDeriv` take vector- and matrix-views, such that we can
take sub-vectors and sub-matrices when calling the evaluations.

We build expression trees, similar as the expression templates for vectors and matrices.
But now we use virtual inheritance instead of the Barton Neckman trick (i.e. dynamic
polymorphism instead of static polymorphism). This is more expensive to create, but it allows to
pass `NonlinearFunction` objects between C++ functions.

A `SumFunction` implements the sum $f_A+f_B$. The two childs are provided by pointers.
Shared pointers allow simple life-time management:

```cpp
class SumFunction : public NonlinearFunction {
  shared_ptr<NonlinearFunction> fa, fb;
public:
  SumFunction (shared_ptr<NonlinearFunction> _fa,
               shared_ptr<NonlinearFunction> _fb)
    : fa(_fa), fb(_fb) { } 
    
  size_t DimX() const { return fa->DimX(); }
  size_t DimF() const { return fa->DimF(); }

  void Evaluate (VectorView<double> x, VectorView<double> f) const
  {
    fa->Evaluate(x, f);
    Vector<> tmp(DimF());
    fb->Evaluate(x, tmp);
    f += tmp;
  }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    fa->EvaluateDeriv(x, df);
    Matrix<> tmp(DimF(), DimX());
    fb->EvaluateDeriv(x, tmp);
    df += tmp;
  }
};
```

To generate such a `SumFunction` object, we overload the `operator+` for two `NonlinearFunction` objects,
represented as shared objectes:
```cpp
auto operator+ (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
{
  return make_shared<SumFunction>(fa, fb);
}
```

## Implementing a Newton solver

Newton's method for solving the non-linear equation

$$
f(x) = 0
$$

is this iterative method:

$$
x^{n+1} = x^n - f^\prime(x^n)^{-1} f(x^n)
$$

If the Jacobi-matrix at the solution $x^\ast$ is regular, and the initial guess $x^0$ is sufficiently close to $x^\ast$,
Newton's method converges quadratically:

$$
\| x^{n+1} - x^\ast \| \leq c \, \| x^n - x^\ast \|^2
$$
This means the number of valid digits double in every iteration.

```cpp
void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                   double tol = 1e-10, int maxsteps = 10,
                   std::function<void(int,double,VectorView<double>)> callback = nullptr)
{
  Vector<> res(func->DimF());
  Matrix<> fprime(func->DimF(), func->DimX());

  for (int i = 0; i < maxsteps; i++)
    {
      func->Evaluate(x, res);
      func->EvaluateDeriv(x, fprime);
      CalcInverse(fprime);
      x -= fprime*res;

      double err = Norm(res);
      if (callback) callback(i, err, x);
      if (err < tol) return;
    }

    throw std::domain_error("Newton did not converge");
  }

```

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