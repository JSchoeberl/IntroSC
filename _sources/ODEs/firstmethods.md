# Some simple time-stepping methods

The goal is to find an approximative solution to the ODE

$$
y^\prime(t) = f(t, y(t)) \; \; \forall \, t \in (t_0, T), \quad y(t_0) = y_0
$$

by a numerical method.

We introduce time points

$$
t_0 < t_1 < \ldots < t_n = T,
$$

usually we choose $n \in {\mathbb N}$, and define $h := \frac{T-t_0}{n}$, and
$t_i := t_0 + ih$.

## Explicit Euler method

We compute a discrete approximative solution 

$$
y_i \approx y(t_i) \qquad \forall \, i \in { 0, \ldots , n}
$$


The derivative $y^\prime$ is approximated by the forward difference quotient

$$
y^\prime (t_i)
\approx \frac{y(t_{i+1}) - y(t_i)}{t_{i+1} - t_i}
\approx \frac{y_{i+1} - y_i}{t_{i+1}-t_i}
$$

Inserting this approximation into the equation we get the numerical scheme

$$
\frac{y_{i+1} - y_i}{t_{i+1}-t_i} = f(t_i, y_i)
$$

or (with constant mesh size $h = t_{i+1} - t_i$)

$$
y_{i+1} = y_i + h f(t_i, y_i) \qquad 0 \leq i < n
$$

This method is very simple, it is a fully explicit algorithm. The next step $y_{i+1}$ is
defined by simply evaluating this formula.

## Implicit Euler method

One can also start with the backward difference formula. The same difference quotient is
now seen as an approximation to $y^\prime(t_{i+1})$. This leads to the discrete scheme

$$
\frac{y_{i+1} - y_i}{t_{i+1}-t_i} = f(t_{i+1}, y_{i+1})
$$

or

$$
y_{i+1} = y_i + h f(t_{i+1}, y_{i+1}) \qquad 0 \leq i < n
$$

Now, the new step $y_{i+1}$ shows up also in the right-hand side. We cannot simply
evaulate a formula, but have to solve a system of equation. If the right hand side $f(t,y)$
depends non-linear on the the state $y$, one has to solve a non-linear system. Typically,
this is done using Newton's method.

## Crank-Nicolson method

Interpreting the difference quotient as an approximation to the derivative in
in the mid-point $\frac{t_{i}+t_{i+1}}{2}$, we end up in the Crank-Nicolson method

$$
y_{i+1} = y_i + \frac{h}{2} \big( f(t_i, y_i) + f(t_{i+1}, y_{i+1}) \big) \qquad 0 \leq i < n
$$

Like the implicit Euler method, it is an implicit method as well. Numerical analysis proves that
this method provides a higher order of convergence.


## Derivation from the integral equation

We consider the ODE

$$
y^\prime(t) = f(t, y(t)) \qquad y(t_i) = y_i
$$

on the time interval $[t_i, t_{i+1}]$. We rewrite it as integral equation

$$
y(t) = y_i + \int_{t_i}^t f(s, y(s)) \, ds \qquad \forall \, t \in [t_i, t_{i+1}]
$$

We insert $t = t_{i+1}$ to obtain the next time-step, 
and use numerial integration to approximate the integral on the right hand side.
The left-sided rectangular rule

$$
\int_a^b f(x) \, dx \approx (b-a) \, f(a)
$$

leads to

$$
y_{i+1} = y_i + h f(t_i, y_i),
$$

the right-sided rectangular rule

$$
\int_a^b f(x) \, dx \approx (b-a) \, f(b)
$$

leads to

$$
y_{i+1} = y_i + h f(t_{i+1}, y_{i+1}).
$$

The more accurate trapezoidal rule

$$
\int_a^b f(x) \, dx \approx (b-a) \, \frac{ f(a) + f(b) }{2} 
$$

leads to the Crank-Nicolson method

$$
y_{i+1} = y_i + \frac{h}{2} \big( f(t_i, y_i) +  f(t_{i+1}, y_{i+1})\big)
$$

We see that we have recovered the same three methods from numerical integration, as we
had already motivated before from numerical differentiation.
We will use more accurate integration rules to develop more accurate time-stepping methods.


