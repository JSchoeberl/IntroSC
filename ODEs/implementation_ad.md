# Automatic Differentiation

Often it is cumbersome and error prune to compute derivatives exactly. Two *different* ways to avoid manual differentiation are numerical differentiation, and automatic differentiation. We don't cover symbolic differentiation as provided by computer algebra libraries as *sympy*.

## Numerical differentiation

To approximate the derivative of a function one may use finite differences, for example use a central difference quotient for the $i^{th}$ partial derivative:

$$
\frac{\partial f}{\partial x_i}(x) \approx \frac{f(x + \varepsilon e_i) - f(x - \varepsilon e_i)}{2 \varepsilon}
$$

Here, $e_i$ is the $i^{th}$ unit-vector, and $\varepsilon$ is a small value. If the function is smooth, one easily verifies using Taylor expansion that the approximation error is $O(\varepsilon^2)$. To have a good approximation, one should choose $\epsilon$ small.

On the other hand, one subtracts two similar values in the numerator, and observes numerical cancellation if $\varepsilon$ is very small. Thus, a good compromise is needed.

The advantage of numerical differentiation is that it is non-invasive, it works without any change of the code evaluating the function $f$.

## Automatic differentiation

With automatic differentiation, we represent the value, and the derivative in every step of the evaluation of a function. From rules of calculus, we know how derivatives are computed for all (differentiable) elementary operations.

As an example, let

$$
f(u,v) = u * v.
$$

Say $(x_1, \ldots x_n)$ are the independent variables after whom we want to differentiate. Say $u$ and $v$ are intermediate results, for which we know the value, and all partial derivatives (the gradient $\nabla u$). They are represented as

$$
{\bf u} = \left( \begin{array}{c} u \atop g_u \end{array} \right) \quad \text{and} \quad
{\bf v} = \left( \begin{array}{c} v \atop g_v \end{array} \right)
$$

Then, by the product rule, we know that

$$
{\bf w} = f( {\bf u}, {\bf v})
$$

should have value and gradient:

$$
{\bf w} = \left( \begin{array}{c} w \atop g_w \end{array} \right) =  \left( \begin{array}{c} u v \atop u g_v + v g_u \end{array} \right)
$$

We now from the chain rule that

$$
\nabla \sin (u) = \cos(u) \nabla u
$$

Thus, the value-gradient tuple is defined as

$$
\sin({\bf u}) = 
\left( \begin{array}{c} \sin(u) \atop \cos(u) g_u \end{array} \right)
$$

The advantage of automatic differentiation is that no choice of a small parameter $\varepsilon$ is needed. However, one needs access to the code implementing the function $f$.

## Implementation of an `AutoDiff` data type

We implement a data type for automatic differentiation. It stores the function value (`m_val`), and the gradient (`m_deriv`). We use a template argument `N` to specify the number of independent variables, what is the dimension of the gradient:

```cpp
template <size_t N, typename T = double>
class AutoDiff
{
  T m_val;
  std::array<T, N> m_deriv;
public: 
  AutoDiff () : m_val(0), m_deriv{} {}
  AutoDiff (T v) : m_val(v), m_deriv{} {}
    
  template <size_t I>
  AutoDiff (Variable<I, T> var) : m_val(var.value()), m_deriv{} {
      m_deriv[I] = 1.0;
  }

  T value() const { return m_val; }
  std::array<T, N>& deriv() { return m_deriv; }
  const std::array<T, N>& deriv() const { return m_deriv; }
};
```

For that class we implement the mathematical functions and operators we want to support, for example

```cpp
template <size_t N, typename T = double>
AutoDiff<N, T> operator* (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b) 
{
   AutoDiff<N, T> result(a.value() * b.value());
   for (size_t i = 0; i < N; i++)
      result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
   return result;
}
```

or

```cpp
template <size_t N, typename T = double>
AutoDiff<N, T> sin(const AutoDiff<N, T> &a)
{
  AutoDiff<N, T> result(sin(a.value()));
  for (size_t i = 0; i < N; i++)
     result.deriv()[i] = cos(a.value()) * a.deriv()[i];
  return result;
}
```

The `Variable<I>` class is used to define and initialize 
the `I`-th independent variable.

It can be used as follows. Say, we have the function

$$
f(x,y) = x * y
$$

and we want to compute

$$
\nabla f(x,y) \quad \text{for} \quad x=1, y=2
$$

We declare our independent variables

```cpp
AutoDiff<2> x = Variable<0>(1.0);
AutoDiff<2> y = Variable<1>(2.0);
```

and define

```cpp
f = x*y;
```

to obtain the gradient.

Find the code in `src\autodiff.hpp` and some demos in `demos\demo_autodiff.cpp`.


## Using numerical and automatic differentiation for the nonlinear-function class

```cpp
class NonlinearFunctionNumDiff : public NonlinearFunction {
  double m_eps;
public:
  NonlinearFunctionNumDiff (double eps) : m_eps(eps) { }
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override {
    Vector xl(dimX()), xr(dimX());
    Vector fl(dimF()), fr(dimF());
    for (int i = 0; i < dimX(); i++) {
      xl = x; xr = x;
      xl(i) -= eps; xr(i) += eps;
      evaluate(xl,fl);
      evaluate(xr,fr);
      df.Col(i) = 1/(2*eps)*(fr-fl);
    }
  };
}
```


```cpp
template <typename NLF>
class NonlinearFunctionAutoDif : public NonlinearFunction {
public:

  void evaluate(VectorView<double> x, VectorView<double> f) const override {
    static_cast<const NFL*>(this) -> T_evaluate(x, f);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override {
    Vector<AutoDiff<1>> adx(dimX());
    Vector<AutoDiff<1>> adf(dimF());

    for (int i = 0; i < dimX(); i++) {
      for (int j = 0; j < dimX(); j++)
        adx(j) = x(j);
      adx(i) = Variable<0>(x(i));
      static_cast<const NFL*>(this) -> T_evaluate(adx, adf);
      for (int j = 0; j < dimF(); j++)
        df(j,i) = adf(j).deriv()[0];
    }
  }  
};
```

and the actual function:

```cpp
class MyFunc : public NonlinearFunctionAutoDiff<MyFunc> {

 public:
 template <typename T>
 void T_evaluate(VectorView<T> x, VectorView<T> f) const {
    // evaluation of function
 }
}
```




## Exercises

* Add additional useful operators for the `AutoDiff` class
* Add some more functions (`cos`, `exp`, `log`, ...) for the `AutoDiff` class.

* Evaluate and plot Legendre-polynomials up to order 5, in the interval $-1 \leq x \leq 1$. Evaluate and plot also
their derivatives (using AutoDiff).

  Legendre polynomials are recursively defined:
  \begin{eqnarray*}
  P_0(x) & = & 1 \\
  P_1(x) & = & x \\
  P_n(x) & = & \frac{2n-1}{n} x P_{n-1}(x) + \frac{n-1}{n} P_{n-2}(x) \qquad n \geq 2
  \end{eqnarray*}

  ```cpp
  template <typename T>
  void LegendrePolynomials(int n, T x, std::vector<T>& P) {
      if (n < 0) {
          P.clear();
          return;
      }
      P.resize(n + 1);
      P[0] = T(1);
      if (n == 0) return;
      P[1] = x;
      for (int k = 2; k <= n; ++k) {
          P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
      }
  }
  ```


## Exercise: test the `AutoDiff` class for the pendulum

The pendulum is modeled by an ode whith non-linear right-hand side. The ode as second order equation is

$$
\alpha^{\prime \prime} = -\frac{g}{l} \sin (\alpha), \qquad \alpha(t_0) = \alpha_0, \; \alpha^\prime(0) = \alpha_{0p},
$$

where $\alpha$ is the angular displacement, $l$ is the length of the pendulum, and $g$ the gravitational acceleration
(roughly $9.81 m s^{-2}$ on earth).

Reformulate the equation as first order ODE in ${\mathbb R}^2$.

We want to implement the right-hand side function only once, where the data type is the template typename `T`,
see function `T_evaluate`. To evaluate $f$, we call the template with a `double` type. To evaluate $f^\prime$,
we call `T_evaluate` with an `AutoDiff<2>` type.


```cpp
class PendulumAD : public NonlinearFunction
{
private:
  double m_length;
  double m_gravity;

public:
  PendulumAD(double length, double gravity=9.81) : m_length(length), m_gravity(gravity) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    T_evaluate<double>(x, f);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    Vector<AutoDiff<2>> x_ad(2);
    Vector<AutoDiff<2>> f_ad(2);

    x_ad(0) = Variable<0>(x(0));
    x_ad(1) = Variable<1>(x(1));
    T_evaluate<AutoDiff<2>>(x_ad, f_ad);

    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
         df(i,j) = f_ad(i).deriv()[j];
  }

  template <typename T>
  void T_evaluate (VectorView<T> x, VectorView<T> f) const

  {
    f(0) = x(1);
    f(1) = -m_gravity/m_length*sin(x(0));
  }
};
```
