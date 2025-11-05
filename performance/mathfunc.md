# Vectorizing mathematical functions

The aim is to implement mathematical functions like `exp`, `sin` and `cos` using the `SIMD` types.


References:

*  Stephen L. Moshier: Methods and Programs For Mathematical Functions
   https://www.moshier.net/methprog.pdf
    
*  CEPHES mathematical function library
   https://www.netlib.org/cephes/


## How to implement `sin` and `cos` ?

* Local approximation by polynomials

  For small values of the argument $x$, polynomial approximation via Taylor expansion can be used.
  A uniform best approxiamtion is obtained via Chebyshev interpolation.
  We took the coefficients from CEPHES, file cmath.tgz, which is the basis of many implementations
     of mathematical functions:
  

  ```cpp
  static constexpr double sincof[] = {
      1.58962301576546568060E-10,
      -2.50507477628578072866E-8,
      2.75573136213857245213E-6,
      -1.98412698295895385996E-4,
      8.33333333332211858878E-3,
      -1.66666666666666307295E-1,
    };
  
  static constexpr double coscof[6] = {
      -1.13585365213876817300E-11,
      2.08757008419747316778E-9,
      -2.75573141792967388112E-7,
      2.48015872888517045348E-5,
      -1.38888888888730564116E-3,
      4.16666666666665929218E-2,
    };
  
  // highly accurate on [-pi/4, pi/4]
  template <typename T>
  auto sincos_reduced (T x) 
  {
    auto x2 = x*x;
    
    auto s = (((((sincof[0]*x2+sincof[1])*x2+sincof[2])*x2+sincof[3])*x2+sincof[4])*x2+sincof[5]);
    s = x + x*x*x * s;
  
    auto c = (((((coscof[0]*x2+coscof[1])*x2+coscof[2])*x2+coscof[3])*x2+coscof[4])*x2+coscof[5]);
    c = 1.0 - 0.5*x2 + x2*x2*c;
  
    return std::tuple{ s, c };
  }
  ```

  This polynomial approximation can be used the same way for `double x`, and `SIMD<double> x`.


* Global approximation

  We represent
  
  $$
  x = q \frac{\pi}{2} + \tilde x
  $$
  with an integer $q$ and the reminder $\tilde x \in [-\pi/4, \pi/4]$.

  Depending on the reminder, we can reduce the evaluation of $\sin x$ to
  local evaluation of $\sin$ or $\cos$ for small $\tilde x$:
  
  $$
  \sin(x) = \left\{ \begin{array}{cc}
          \sin(\tilde x) & \text{for } q \text{ mod } 4 = 0, \\
          \cos(\tilde x) & \text{for } q \text{ mod } 4 = 1, \\
          -\sin(\tilde x) & \text{for } q \text{ mod } 4 = 2, \\
          -\cos(\tilde x) & \text{for } q \text{ mod } 4 = 3 
          \end{array} \right.
  $$
  
  The choice of the branch can be made by the conditional operator `?` applied to the lowest bits of `q`:
  
  ```cpp
  auto sincos (double x)
  {
    double y = round((2/M_PI) * x);
    int q = lround(y);
    
    auto [s1,c1] = sincos_reduced(x - y * (M_PI/2));
  
    double s2 = ((q & 1) == 0) ? s1 : c1;
    double s  = ((q & 2) == 0) ? s2 : -s2;
  
    double c2 = ((q & 1) == 0) ? c1 : -s1;
    double c  = ((q & 2) == 0) ? c2 : -c2;
    
    return std::tuple{ s, c };
  }
  ```


  All these operations can be perfectly simded:
  
  ```cpp
  template <int N>
  auto sincos (SIMD<double,N> x)
  {
    SIMD<double,N> y = round((2/M_PI) * x);
    SIMD<int64_t,N> q = lround(y);
    
    auto [s1,c1] = sincos_reduced(x - y * (M_PI/2)); 
  
    auto s2 = select((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), s1,  c1);
    auto s  = select((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), s2, -s2);
    
    auto c2 = select((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), c1, -s1);
    auto c  = select((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), c2, -c2);
    
    return std::tuple{ s, c };
  }
  ```



## What the f..k?

An important operation in compute graphics (games!) is the normalization of a vector, $|v| = 1/\sqrt(v*v) v$.
This requires the inverse square root function

$$
y = \frac{1}{\sqrt{x}}
$$

Instead of calling a `sqrt` function and a division, one may solve the nonlinear equation

$$
F(y) = \frac{1}{y^2} - x = 0
$$

using a few steps of Newton's method:

\begin{eqnarray*}
y_{n+1} & = & y_n - \frac{F(y_n)}{F^\prime(y_n)} \\
   & = & y_n - \frac{ y_n^{-2} - x } { - 2 y_n^{-3} } \\
   & = & \tfrac{3}{2} y - \tfrac{1}{2} x y_n^2
\end{eqnarray*}

Every step requires only 3 multiplications.

This algorithm got famous in [Quake III Arena](https://en.wikipedia.org/wiki/Fast_inverse_square_root)
providing fluent 3D graphics (in the year 1999).
Find the full code including the documented initial value for Newton's method at this page.

For example, processors by Arm provide intrinsics exactly for that Newton iteration:
[instruction](https://docs.unity3d.com/Packages/com.unity.burst@1.6/api/Unity.Burst.Intrinsics.Arm.Neon.vrsqrtsq_f64.html)

## Exercise

* add missing SIMD-functions for the $\sin$ and $\cos$ function.

  Check accuary of the function by comparing to the standard $\sin$ and $\cos$ functions.
  
  Measure speed by calling $10^8$ functions.


* Implement the $\exp$ function for `SIMD` types.

  Rewrite for 

  $$
  x = q \ln 2 + \tilde x
  $$

  with $q \in {\mathbb Z}$ and reminder $\tilde x \in [-\tfrac{\ln 2}{2}, \tfrac{\ln 2}{2}]$:

  \begin{eqnarray*}
  e^x & = & e^{ q \ln 2 + \tilde x} = e^{q \ln2 } e^\tilde x \\
  & = & 2^q \, e^{\tilde x}
  \end{eqnarray*}

  * implement an approximation of $\exp(x)$ for small $x$

  * to implement $2^q$, use the representation of [floating point numbers](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)


Some places to find insparation:
* [Agner Fogs' vectorclass](https://github.com/vectorclass/version2/blob/master/vectormath_trig.h)
* [SLEEF library](https://sleef.org)
* [Netgen/NGSolve simd functions](https://github.com/NGSolve/netgen/blob/master/libsrc/core/simd_math.hpp)






