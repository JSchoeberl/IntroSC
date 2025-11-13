# Function algebra

We implement a function algebra, which allows us to write expressions like
```cpp
  func = Id + 3 * Compose (f, g);
```
where `f` and `g` are functions, and `Id` is the identic function. Then
the composed function `func` shall be able to compute the function value
and the derivative at a given point $x$ (of type `Vector`).
```cpp
Vector vecf(func->dimF());
func->Evaluate(x, vecf);
 
Matrix jacobi(func->dimF(), func->dimX());
func->EvaluateDeriv(x, jacobi);
```

Our base class for such functions is
```cpp
class NonlinearFunction {
public:
  virtual ~NonlinearFunction() = default;
  virtual size_t dimX() const = 0;
  virtual size_t dimF() const = 0;
  virtual void evaluate (VectorView<double> x, VectorView<double> f) const = 0;
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
};
```
`dimX` and `dimF` provide the vector space dimensions of the domain, and the image.
The `evaluate` and `evaluateDeriv` take vector- and matrix-views, such that we can
take sub-vectors and sub-matrices when calling the evaluations.

The function algebra is implemented in the file `nonlinfunc.hpp`.


We build expression trees, similar as the expression templates for vectors and matrices.
But now we use virtual inheritance instead of the Barton Neckman trick (i.e. dynamic
polymorphism instead of static polymorphism). This is more expensive to create, but it allows to
pass `NonlinearFunction` objects between C++ functions.

A `SumFunction` implements the sum $f_A+f_B$. The two childs are provided by pointers.
Shared pointers allow simple life-time management:

```cpp
class SumFunction : public NonlinearFunction {
  shared_ptr<NonlinearFunction> m_fa, m_fb;
public:
  SumFunction (shared_ptr<NonlinearFunction> fa,
               shared_ptr<NonlinearFunction> fb)
    : m_fa(fa), m_fb(fb) { } 
    
  size_t dimX() const { return m_fa->dimX(); }
  size_t dimF() const { return m_fa->dimF(); }

  void evaluate (VectorView<double> x, VectorView<double> f) const
  {
    m_fa->evaluate(x, f);
    Vector<double> tmp(dimF());
    m_fb->Evaluate(x, tmp);
    f += tmp;
  }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    fa->dvaluateDeriv(x, df);
    Matrix<double> tmp(dimF(), dimX());
    fb->dvaluateDeriv(x, tmp);
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


Some more building blocks for our function hierarchy are:

* The IdentityFunction $f(x) = x$
* A constant function $f(x) = val$
* A scaled function $\alpha f$, wehere $\alpha$ is a scalar
* A composed function $f \circ g$
* ...

