# Runge-Kutta methods

The goal of Runge-Kutta methods is to obtain higher accuarcy by
using a more accurate integration rule:

$$
y_{i+1} = y_i + h \sum_{j=1}^s b_j f(t_i+c_j h, y_i^j)
$$

where

* $s$ is called the number of stages

* $c_j, b_j$ are points and weights of the quadrature rule

$$
\int_0^1 f(t) dt \approx \sum_{j=1}^s b_j f(c_j)
$$

* $y_i^j$ are approximations to $y(t_i + c_j h)$, which are also unknown

The exact values at the stage-points are

$$
y(t_i+c_j h) = y(t_i) + \int_{t_i}^{t_i+c_j h} f(t, y(t)) dt
$$

Again, we use numerical integration to define the numerical stage values (skipping the index $i$):

$$
y^j = y_i + h \sum_{k=1}^s  a_{jk} f(t_i+c_k h, y^k) \qquad 1 \leq j \leq s
$$

We use the same integration points as above, however the weigths are adjusted
to the smaller interval such that

$$
\int_0^{c_j} f(t) dt \approx \sum_{k=1}^s a_{jk} f(c_k) \qquad 1 \leq j \leq s
$$

A particular RK method is completely defined by providing the $a$, $b$, and $c$ coefficients,
which are usually given in the so called **Butcher tableau**:


$$
\begin{array}{c|ccc}
c_1 & a_{11} & \dots & a_{1s} \\
\vdots & \vdots & & \vdots \\
c_s & a_{s1} & \dots & a_{ss} \\
\hline
 & b_1 & \dots & b_s
\end{array}
$$


## Explicit RK methods

Explicit RK methods can be seen as extension of the explicit Euler method.
The $a_{jk}$ coefficients form a strictly lower triangular matrix.

Examples are

* the $RK_2$ - method, aka the explicit mid-point rule

$$
\begin{array}{c|cc}
0 & 0 & 0 \\
\frac{1}{2} & \frac{1}{2} & 0 \\
\hline
 & 0 & \frac{1}{2}
\end{array}
$$


* the $RK_4$ - method, aka THE classical Runge Kutta method

$$
\begin{array}{c|cccc}
0 &   \\
\frac{1}{2} & \frac{1}{2}  \\
\frac{1}{2} & 0 & \frac{1}{2} \\
1 & 0 & 0 & 1 \\
\hline
 & \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6}
\end{array}
$$



## Implicit Runge Kutta methods

Gauss-Legendre 2 method:

$$
\begin{array}{c|cccc}
\frac{1}{2} - \frac{\sqrt{3}}{6} &  \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{6} \\
\frac{1}{2} + \frac{\sqrt{3}}{6} & \frac{1}{4} + \frac{\sqrt{3}}{6}  & \frac{1}{4} \\
\hline
 & \frac{1}{2} & \frac{1}{2}
\end{array}
$$


Radau IIA method

$$
\begin{array}{c|cc}
\frac{1}{3} & \frac{5}{12} & \frac{-1}{12} \\
1 & \frac{3}{4} & \frac{1}{4} \\
\hline
  & \frac{3}{4} & \frac{1}{4} \\
\end{array}
$$


**Exercise**: Implement a Runge-Kutta time-stepper for arbitrary Butcher tableaus.

Implement a class like 
```cpp
class ButcherFunction : public NonlinearFunction {
  Matrix a;  // a-coef of Butcher tableau
  shared_ptr<NonlinearFunction> func;
public:
  virtual Evaluate (VectorView x, VectorView f) {
    f = 0.0;
    for (int j = 0; j < s; j++) {
        func->Evaluate (x.Range(j*n, (j+1)*n), tmp);
        for (int k = 0; k < s; k++)
            f.Range(k*n, (k+1)*n) += a(k,j) * tmp;
  }
}
```
Setup the equation, and solve it using `NewtonSolver`.


### Runge-Kutta methods of arbitrary order

Step 1: Start from an Gaussian numerical integration rule to define $c$ and $b$ coefficients.
Take care to scale to the  interval $[0,1]$.

A ready to use function can be taken from 
[Numerical recipes, section 4.6](https://numerical.recipes/book.html) 


Step 2: Find coefficients of the matrix $a$. Each row should provide an integraton rule
on $[0, c_j]$, being exact for polynomials up to order $s-1$:

$$
\sum_{k=1}^s c_k a_{jk} x^{j-1} = \int_0^{c+j} x^{j-1} \, dx \qquad 1 \leq j \leq s
$$

This gives $s$ conditions for each row. Setup linear systems of equations, and solve them.




