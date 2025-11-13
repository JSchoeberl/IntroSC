# Automatic Differentiation

Often it is cumbersome and error prune to compute derivatives exactly. Two *different* ways to avoid manual differentiation are numerical differentiation, and automatic differentiation.

## Numerical differentiation

To approximate the derivative of a function one may use finite differences, for example use a central difference quotient for the $i^{th}$ partial derivative:

$$
\frac{\partial f}{\partial x_i}(x) \approx \frac{f(x + \varepsilon e_i) - f(x - \varepsilon e_i)}{2 \varepsilon}
$$

Here, $e_i$ is the $i^{th}$ unit-vector, and $\varepsilon$ is a small value. If the function is smooth, one easily verifies using Taylor expansion that the approximation is $O(\varepsilon^2)$. To have a good approximation, one should choose $\epsilon$ small.

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
{\bf w} = \left( \begin{array}{c} w \atop g_u \end{array} \right) =  \left( \begin{array}{c} u v \atop u g_v + v g_u \end{array} \right)
$$

We now that from the chain rule that

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

```cppp
template <size_t N, typename T = double>
class AutoDiff
{
  T m_val;
  std::array<T, N> m_deriv;
public: 
  AutoDiff () : m_val(0), m_deriv{} {}
  AutoDiff (T v) : m_val(v), m_deriv{} 
    
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
the `I`th independent variable.

It can be used as follows. Say, we have the function

$$
f(x,y) = x * y
$$

and we want to compute

$$
\nabla f(x,y) \quad \text{for} \quad x=1, y=2
$$

For that we declare our independent variables

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


## Exercises

* add additional useful operators for the `AutoDiff` class
* add some more functions (`cos`, `exp`, `log`, ...) for the `AutoDiff` class.
