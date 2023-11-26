# Implementation

We implement a function algebra, which allows us to write expressions like
```cpp
  func = Id + 3 * Compose (f, g);
```
where `f` and `g` are functions, and `Id` is the identic function. Then
the composed function `func` shall be able to provide the function value and the
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
Newton's method converges quadratic:

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

## Coding the Emplicit Euler method

In every time-step we have to solve for the new value $y^{n+1}$:

$$
y^{n+1} - y^n - h f(y^{n+1}) = 0
$$

The function $f : {\mathbb R}^n \rightarrow {\mathbb R}^n$ is the right hand side of the ODE,
which has been brought into autonomous form.

We use our function algebra to build this composed function, and throw it into the Newton solver.
If we make the time-step not too large, the value $y^n$ of the old time-step is a good starting value.

To express that the independent varible is $y^{n+1}$, we create an `IdenticFunction`. The old
value is considered as a constant, i.e. a `ConstantFunction`. The right hand side function is
given by the user. Then the implicit Euler method is coded up like that:

```cpp
void SolveODE_IE(double tend, int steps,
                 VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                 std::function<void(double,VectorView<double>)> callback = nullptr)
{
  double dt = tend/steps;
  auto yold = make_shared<ConstantFunction>(y);
  auto ynew = make_shared<IdenticFunction>(y.Size());
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

