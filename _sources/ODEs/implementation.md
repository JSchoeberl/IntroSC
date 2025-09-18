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
