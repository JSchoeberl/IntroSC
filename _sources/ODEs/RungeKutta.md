# Runge-Kutta methods

The goal of Runge-Kutta methods is to obtain higher accuarcy by
using a more accurate integration rule:

$$
y_{i+1} = y_i + h \sum_{j=0}^{s-1} b_j f(t_i+c_j h, y_i^j)
$$

where

* $s$ is called the number of stages

* $c_j, b_j$ are points and weights of the quadrature rule

$$
\int_0^1 f(t) dt \approx \sum_{j=0}^{s-1} b_j f(c_j)
$$

* $y_i^j$ are approximations to $y(t_i + c_j h)$, which are also unknown

The exact values at the stage-points are

$$
y(t_i+c_j h) = y(t_i) + \int_{t_i}^{t_i+c_j h} f(t, y(t)) dt
$$

Again, we use numerical integration to define the numerical stage values (skipping the index $i$):

$$
y^j = y_i + h \sum_{l=0}^{s-1}  a_{jl} f(t_i+c_l h, y^l) \qquad 0 \leq j < s
$$

We use the same integration points as above, however the weigths are adjusted
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

Instead of having the $y^j$ as unknonws, one often solves for $k$-values ($k$ like slope):

$$
k^j = f(y_i + c_j h, y_i + h \sum_{l=0}^{s-1} a_{jl} k^l)
$$

and then

$$
y_{i+1} = y_i + h \sum_{l=0}^{s-1} b_l k^l
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
y^0 & = & y_i \\
y^1 & = & y_i + \tfrac{h}{2} f(t_i, y^0)
\end{array}
$$

and then

$$
y_{i+1} = y_i + h f(t_i+\tfrac{h}{2}, y^1)
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
y^0  & = & y_i \\
y^1  & = & y_i + \tfrac{h}{2} f(t_i, y^0) \\
y^2  & = & y_i + \tfrac{h}{2} f(t_i+\tfrac{h}{2}, y^1) \\
y^3  & = & y_i + h f(t_i+\tfrac{h}{2}, y^2) \\
\end{array}
$$

and then

$$
y_{i+1} = y_i + \frac{h}{6} \big(f(t_i, y^0) + 2 f(t_i+\tfrac{h}{2}, y^1)
+ 2 f(t_i+\tfrac{h}{2}, y^2) + f(t_i+h, y^3) \big)
$$




## Implicit Runge Kutta methods

Methods of optimal accuracy are obtained based on Gaussian quadrature:

Gauss-Legendre 2 method:

$$
\begin{array}{c|cccc}
\frac{1}{2} - \frac{\sqrt{3}}{6} &  \frac{1}{4} & \frac{1}{4} - \frac{\sqrt{3}}{6} \\
\frac{1}{2} + \frac{\sqrt{3}}{6} & \frac{1}{4} + \frac{\sqrt{3}}{6}  & \frac{1}{4} \\
\hline
 & \frac{1}{2} & \frac{1}{2}
\end{array}
$$


If not all time-scales are resolved, then $c_s = 1$ leads to better stability.
The other points are chosen for optimal accuracy. These methods are called Radau IIA.

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

## Exercise: Implement a Runge-Kutta time-stepper for arbitrary Butcher tableaus.

Solve the non-linear system of equation for $k \in {\mathbb R}^{sn}$:

$$
k - \tilde f ( \tilde y_i  + h A \otimes k) = 0
$$

Where
* $k = (k_0, \ldots k_{s-1})$, with $k_j \in {\mathbb R}^n$
* $\tilde y_i = (\underbrace{y_i, \ldots , y_i}_{s-{\text times}}) \in {\mathbb R}^{sn}$
* $\tilde f : {\mathbb R}^{sn} \rightarrow {\mathbb R}^{sn} : (x_0, \ldots, x_{s_-1}) \mapsto (f(x_0), \ldots, f(x_{s-1}))$
* for $A \in {\mathbb R}^{h \times w}$ define
$A \otimes :  {\mathbb R}^{wn} \rightarrow {\mathbb R}^{hn} : (x_0, \ldots, x_{w-1}) \mapsto (y_0, \ldots , y_{h-1})$ such that $y_j = \sum_l a_{jl} x_l$.

Derive two classes from `NonlinearFunction` for that (for example, named `BlockFunction` and `BlockMatVec`).

Setup the Runge-Kutta equation, and solve it using `NewtonSolver`.


### Runge-Kutta methods of arbitrary order

If we have fixed the integration points $c_j$, we can compute the $a_{jl}$ and $b_j$ coefficients
such that the integration rules are (at least) exact for polynomials up to order $s-1$:

$$
\sum_{l=0}^{s-1} b_{l} {c_l}^k = \int_0^1 x^k \, dx \qquad 0 \leq k < s
$$

and

$$
\sum_{l=0}^{s-1} a_{jl} {c_l}^k = \int_0^{c_j} x^k \, dx \qquad 0 \leq j,k < s
$$

Setup and solve linear equations for computing the $a$ and $b$ coefficients.


Use the $c$-coefficients from Gauss-Legendre 2, and two-point Radau IIA, and check whether you recover the same $a$ and $b$ coefficients.


Try $c$-coefficients from Gauss-Legendre 3, and three-point Radau II A (from the web-page above).


Then, find $c$-coefficients by solving for the roots of Legendre-polynomials, and transform to $[0,1]$.
A ready to use function can be taken from 
[Numerical recipes, section 4.6](https://numerical.recipes/book.html) 







