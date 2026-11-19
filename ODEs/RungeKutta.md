# Runge-Kutta methods

Recall the integral formulation of the ODE on time interval $[t_n, t_{n+1}]$:

$$
y(t_{n+1}) = y(t_n) + \int_{t_n}^{t_{n+1}} f(s, y(s)) ds
$$

The goal of Runge-Kutta methods is to obtain higher accuracy by
using more accurate numerical integration rule:

$$
y_{n+1} = y_n + \tau \sum_{j=0}^{s-1} b_j f(t_n+c_j \tau, y_n^j)
$$

where

* $s$ is called the number of stages

* $c_j, b_j$ are points and weights of the quadrature rule

$$
\int_0^1 f(t) dt \approx \sum_{j=0}^{s-1} b_j f(c_j)
$$

* $y_n^j$ are approximations to $y(t_n + c_j \tau)$, which are also unknown

The exact values at the stage-points are

$$
y(t_n+c_j \tau) = y(t_n) + \int_{t_n}^{t_n+c_j \tau} f(t, y(t)) dt
$$

Again, we use numerical integration to define the numerical stage values (skipping the index $n$):

$$
y^j = y_n + \tau \sum_{l=0}^{s-1}  a_{jl} f(t_n+c_l \tau, y^l) \qquad 0 \leq j < s
$$

We use the same integration points as above, however the weights are adjusted
to the smaller interval such that

$$
\int_0^{c_j} f(t) dt \approx \sum_{l=0}^{s-1} a_{jl} f(c_l) \qquad 0 \leq j < s
$$

A particular RK method is completely defined by providing the $a$, $b$, and $c$ coefficients,
which are usually given in the so called **Butcher tableau**:


$$
\begin{array}{c|ccc}
c_0 & a_{0,0} & \dots & a_{0,s-1} \\
\vdots & \vdots & & \vdots \\
c_{s-1} & a_{s-1,0} & \dots & a_{s-1,s-1} \\
\hline
 & b_0 & \dots & b_{s-1}
\end{array}
$$

Instead of having the $y^j$ as unknowns, one often solves for $k$-values ($k$ like slope):

$$
k^j = f(t_n + c_j \tau, y_n + \tau \sum_{l=0}^{s-1} a_{jl} k^l)
$$

and then

$$
y_{n+1} = y_n + \tau \sum_{l=0}^{s-1} b_l k^l
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
   & 0 & 1
  \end{array}
  $$

  It is evaluated as

  $$
  \begin{array}{rcl}
  y^0 & = & y_n \\
  y^1 & = & y_n + \tfrac{\tau}{2} f(t_n, y^0)
  \end{array}
  $$

  and then

  $$
  y_{n+1} = y_n + \tau f(t_n+\tfrac{\tau}{2}, y^1)
  $$

  Evaluation in $k$-values:

  $$
  \begin{array}{rcl}
  k^0 & = & f(t_n, y_n) \\
  k^1 & = & f(t_n+\frac{\tau}{2}, y_n + \frac{\tau}{2} k^0)
  \end{array}
  $$
  
  and then

  $$
  y_{n+1} = y_n + \tau k^1
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


$$
\begin{array}{rcl}
y^0  & = & y_n \\
y^1  & = & y_n + \tfrac{\tau}{2} f(t_n, y^0) \\
y^2  & = & y_n + \tfrac{\tau}{2} f(t_n+\tfrac{\tau}{2}, y^1) \\
y^3  & = & y_n + \tau f(t_n+\tfrac{\tau}{2}, y^2) \\
\end{array}
$$

and then

$$
y_{n+1} = y_n + \frac{\tau}{6} \big(f(t_n, y^0) + 2 f(t_n+\tfrac{\tau}{2}, y^1)
+ 2 f(t_n+\tfrac{\tau}{2}, y^2) + f(t_n+\tau, y^3) \big)
$$

  
 
  $$
  \begin{array}{rcl}
  k^0  & = & f(t_n, y_n) \\
  k^1  & = & f(t_n+\frac{\tau}{2}, y_n + \frac{\tau}{2} k^0) \\
  k^2  & = & f(t_n+\frac{\tau}{2}, y_n + \frac{\tau}{2} k^1) \\
  k^3  & = & f(t_n + \tau, y_n + \tau k^2)
  \end{array}
  $$

  and then

  $$
  y_{n+1} = y_n + \frac{\tau}{6} \big( k^0 + 2 k^1 + 2 k^2 + k^3 \big)
  $$






## Implicit Runge Kutta methods

Methods of optimal accuracy (order $2s$) are obtained based on Gaussian quadrature:

Gauss-Legendre 2 method:

$$
\begin{array}{c|cc}
\frac{1}{2} - \frac{\sqrt{3}}{6} &  \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{6} \\
\frac{1}{2} + \frac{\sqrt{3}}{6} & \frac{1}{4} + \frac{\sqrt{3}}{6}  & \frac{1}{4} \\
\hline
 & \frac{1}{2} & \frac{1}{2}
\end{array}
$$


If not all time-scales are resolved, then $c_{s-1} = 1$ leads to better stability (i.e. L-stability).
The other points are chosen for optimal accuracy. These methods are called Radau IIA, they have order $2s-1$.

Two-point Radau IIA method:

$$
\begin{array}{c|cc}
\frac{1}{3} & \frac{5}{12} & \frac{-1}{12} \\
1 & \frac{3}{4} & \frac{1}{4} \\
\hline
  & \frac{3}{4} & \frac{1}{4} \\
\end{array}
$$

Many examples can be found here: [Butcher tableaus examples](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods).


(sec:rk-order)=
## Order and stability of Runge-Kutta methods

We always assume the row-sum condition $c_j = \sum_l a_{jl}$, i.e. the stage $j$ integrates
exactly from $t_n$ to $t_n + c_j \tau$ for constant functions.
By Taylor expansion of the exact and the numerical solution one derives conditions
on the coefficients such that the local error is of order $\tau^{p+1}$, i.e. the method is
consistent of order $p$ (see section {ref}`sec:convergence`):

| order | conditions |
|-------|------------|
| $p \geq 1$ | $\sum_j b_j = 1$ |
| $p \geq 2$ | $\sum_j b_j c_j = \frac{1}{2}$ |
| $p \geq 3$ | $\sum_j b_j c_j^2 = \frac{1}{3}, \quad \sum_{j,l} b_j a_{jl} c_l = \frac{1}{6}$ |
| $p \geq 4$ | $\sum_j b_j c_j^3 = \frac{1}{4}, \quad \sum_{j,l} b_j c_j a_{jl} c_l = \frac{1}{8}, \quad \sum_{j,l} b_j a_{jl} c_l^2 = \frac{1}{12}, \quad \sum_{j,l,m} b_j a_{jl} a_{lm} c_m = \frac{1}{24}$ |

The first condition in every row says that the quadrature rule $(b, c)$ is exact for polynomials of degree $p-1$.
The other conditions involve the $a$-coefficients, they make sure that the stage values are accurate enough.
The number of conditions grows quickly with the order. An explicit method with $s$ stages has order
at most $s$, and for order $p \geq 5$ one needs more than $p$ stages.

Applied to the test equation $y^\prime = \lambda y$, the stage equations become linear, and we obtain
the stability function (see section {ref}`sec:stability`)

$$
g(z) = 1 + z \, b^T (I - z A)^{-1} {\mathbb 1}, \qquad {\mathbb 1} = (1, \ldots, 1)^T.
$$

For explicit methods $A$ is strictly lower triangular, and $g$ is a polynomial. Thus $|g(z)| \rightarrow \infty$ for $|z| \rightarrow \infty$,
and no explicit method is $A$-stable. For $s = p \leq 4$ the stability function is the Taylor polynomial
$g(z) = \sum_{k=0}^{s} z^k / k!$ of $e^z$.

| method | stages $s$ | order $p$ | type | $A$-stable | $L$-stable |
|--------|-----------|-----------|------|------------|------------|
| explicit Euler | 1 | 1 | explicit | no | no |
| explicit mid-point ($RK_2$) | 2 | 2 | explicit | no | no |
| classical $RK_4$ | 4 | 4 | explicit | no | no |
| implicit Euler (= Radau IIA, $s=1$) | 1 | 1 | implicit | yes | yes |
| implicit mid-point (= Gauss-Legendre, $s=1$) | 1 | 2 | implicit | yes | no |
| Crank-Nicolson | 2 | 2 | implicit | yes | no |
| Gauss-Legendre | $s$ | $2s$ | implicit | yes | no |
| Radau IIA | $s$ | $2s-1$ | implicit | yes | yes |

The Crank-Nicolson method is the RK method with $c = (0, 1)$, $b = (\frac{1}{2}, \frac{1}{2})$, and $a_{0,\cdot} = (0, 0)$, $a_{1,\cdot} = (\frac{1}{2}, \frac{1}{2})$.

## Implementation of a Runge-Kutta time-stepper for arbitrary Butcher tableaus.

They are implemented in ASC-ODE, branch RungeKutta. Merge that branch into your main branch:

    git remote add upstream https://github.com/TUWien-ASC/ASC-ODE.git
    git fetch upstream
    git merge upstream/RungeKutta


For an autonomous ODE, we solve a non-linear system of equations for $k^0, \ldots k^{s-1}$ such that

$$
\left( \begin{array}{c}
k^0 \\ \vdots \\ k^{s-1}
\end{array} \right)
-
\left( \begin{array}{c}
f( y_n + \tau (a_{0,0} k^0 + \cdots + a_{0,s-1} k^{s-1}) ) \\
\vdots \\
f( y_n + \tau (a_{s-1,0} k^0 + \cdots + a_{s-1,s-1} k^{s-1}) )
\end{array} \right)
= 0
$$

This system of equations is rewritten in compact form for
the unknown $k \in {\mathbb R}^{sd}$:

$$
k - \widetilde f ( \tilde y_n  + \tau A \otimes k) = 0
$$

Where
* $k = (k^0, \ldots k^{s-1})$, with $k^j \in {\mathbb R}^d$
* $\tilde y_n = (\underbrace{y_n, \ldots , y_n}_{s \text{ times}}) \in {\mathbb R}^{sd}$
* $\tilde f : {\mathbb R}^{sd} \rightarrow {\mathbb R}^{sd} : (x_0, \ldots, x_{s-1}) \mapsto (f(x_0), \ldots, f(x_{s-1}))$
* for $A \in {\mathbb R}^{p \times q}$ define
$A \otimes :  {\mathbb R}^{qd} \rightarrow {\mathbb R}^{pd} : (x_0, \ldots, x_{q-1}) \mapsto (z_0, \ldots , z_{p-1})$ such that $z_j = \sum_l a_{jl} x_l$,
i.e. multiplication with the Kronecker product $A \otimes I_d$.

We derive two classes from `NonlinearFunction` for that (named `MultipleFunc` and `MatVecFunc`).
With them, the class `ImplicitRungeKutta` in `src/implicitRK.hpp` sets up the equation for $k$ by the function algebra:

```cpp
m_equ = knew - Compose(multiple_rhs, m_yold + m_tau * std::make_shared<MatVecFunc>(a, m_n));
```

and solves it with Newton's method in every time-step.



### Runge-Kutta methods of arbitrary order

If we have fixed the integration points $c_j$, we can compute the $a_{jl}$ and $b_j$ coefficients
such that the integration rules are (at least) exact for polynomials up to degree $s-1$:

$$
\sum_{l=0}^{s-1} b_{l} {c_l}^k = \int_0^1 x^k \, dx \qquad 0 \leq k < s
$$

and

$$
\sum_{l=0}^{s-1} a_{jl} {c_l}^k = \int_0^{c_j} x^k \, dx \qquad 0 \leq j,k < s
$$

These are linear equations with the Vandermonde matrix $(c_l^k)_{k,l}$, they are solved by the function `computeABfromC`.
The resulting methods are called collocation methods: the stage values are the values of the polynomial of degree $s$
which satisfies the ODE in the points $t_n + c_j \tau$.

If the $c_j$ are the Gauss points, the quadrature rule $(b,c)$ is exact even up to degree $2s-1$,
and we obtain the Gauss-Legendre methods of order $2s$.
The Gauss points are the roots of the Legendre polynomial $P_s$, transformed from $[-1,1]$ to $[0,1]$.

For the Radau IIA methods, the points are the roots of the Jacobi polynomial $P^{(1,0)}_{s-1}$, transformed to $[0,1]$,
together with $c_{s-1} = 1$.

The functions `GaussLegendre` and `GaussRadau` compute these points, they are taken from
[Numerical recipes, section 4.6](https://numerical.recipes/book.html).


## Exercises

* Try out and compare various Runge-Kutta methods already implemented for the examples above.

* Use `computeABfromC` with the $c$-coefficients from Gauss-Legendre 2 and two-point Radau IIA, and check whether you recover the
  Butcher tableaus above. Do the same for Gauss-Legendre 3 and three-point Radau IIA (from the web-page above).

* Measure the order of convergence of Gauss-Legendre and Radau IIA methods with $s = 1, 2, 3$ stages, as in
  the exercise for the implicit Euler method. Do you observe orders $2s$ and $2s-1$?
  What happens for very small time-steps?

* Compute the Gauss points yourself: evaluate $P_s$ and its derivative with the Legendre polynomials
  and the `AutoDiff` class from section {doc}`implementation_ad`, and find the roots by Newton's method.
  Good initial guesses are $x_j = \cos \big( \pi \frac{j + 3/4}{s + 1/2} \big)$, $0 \leq j < s$. Compare with `GaussLegendre`.

* Implement and test an `ExplicitRungeKutta` time-stepper for arbitrary Butcher tableaus (with a strictly lower triangular matrix for
the $a$-coefficients). Verify the order conditions for the $RK_4$ method.


