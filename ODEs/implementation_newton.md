# Implementing a Newton solver

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

Here is a simple implementation of Newton's method using `NonlinearFunction`.
The vector `x` is used to provide the initial guess, as well as for the solution:

```cpp
void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                   double tol = 1e-10, int maxsteps = 10,
                   std::function<void(int,double,VectorView<double>)> callback = nullptr)
{
  Vector<> res(func->dimF());
  Matrix<> fprime(func->dimF(), func->dimX());

  for (int i = 0; i < maxsteps; i++)
    {
      func->evaluate(x, res);
      func->evaluateDeriv(x, fprime);
      calcInverse(fprime);
      x -= fprime*res;

      double err = norm(res);
      if (callback) callback(i, err, x);
      if (err < tol) return;
    }

    throw std::domain_error("Newton did not converge");
  }

```


It can be used as follows:
```cpp
shared_ptr<NonlinearFunction> myfunc = ...
Vector<double> x = { 1.0, 2.0, 3.0 };

NewtonSolver (myfunc, x, 1e-8, 20,
             [](int step, double err, VectorView<double> x) {
             cout << "step " << step << ", err " << err << ", x = " << x << endl;
             }
```

