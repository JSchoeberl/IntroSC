
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


