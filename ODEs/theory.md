# A little bit of theory

An ordinary differential equation (ODE) is given by

$$
y^\prime (t) = f(t,y(t)) \qquad \forall \, t \in (t_0, T)
$$
together with the initial condition

$$
y(t_0) = y_0,
$$

where
* $f : [t_0, T] \times {\mathbb R}^n \rightarrow {\mathbb R}^n$ is the given **right hand side**
* $y : [t_t, T] \rightarrow {\mathbb R}^n$ is the unknown function, called **state**
* $t \in [t_0, T]$ is called **time**
* $y_0 \in {\mathbb R}^n$ is the **initial value**


## Some examples:

* $n = 1$, $f(t,y) = a y$:

$$
y^\prime(t) = a y, \quad y(0) = y_0
$$

This ode has the solution $y(t) = y_0 e^{a t}$. If $a > 0$, the solution is exponentiall increasing,
for $a < 0$, it is exponentially decreasing.


* Consider $y^{\prime \prime}(t) = -\omega^2 y(t)$ with $y(0) = y_0$, $y^\prime(0) = y_{0p}$. This is a second order
ODE, which can be reduced to a first oder system with the definition $y_1 := y$ and $y_2 := y^\prime$:

$$
\left( \begin{array}{c} y_1 \\ y_2 \end{array} \right)^\prime =
\left( \begin{array}{c} y_2 \\ -\omega^2 y_1 \end{array} \right) 
\quad \text{with i.c.} \quad
\left( \begin{array}{c} y_1 \\ y_2 \end{array} \right)(0) =
\left( \begin{array}{c} y_0 \\ y_{0p} \end{array} \right) 
$$

Its solution is $y(t) = y_0 \cos( \omega t ) + \frac{y_{0p}}{\omega} \sin (\omega t)$.


## Some Theorems

**Theorem: Picard Lindelöf**

Assume the right hand side satisfies a Liptschitz condition in the second argument

$$
\| f(t, y) - f(t, z) \| \leq L \, \| y - z \| \qquad \forall \, t \in [t_0, T], \, \forall \, y,z \in {\mathbb R}^n.
$$

Then there exisits a unique solution of the ODE.

**Theorem: Stability with respect to initial conditions**

Consider two ODEs differing only in the initial values:

$$
y^\prime(t) = f(t, y(t)), \quad y(t_0) = y_0 \\
$$
$$
z^\prime(t) = f(t, z(t)), & & \quad z(t_0) = z_0 \\
$$

Then there holds:

$$
\| y(t) - z(t) \| \leq e^{L (t-t_0)} \| y_0 - z_0 \|
$$

This means that an error in the initial condition may be propagated with exponential groth with constant $L$.


**Theorem: One-sided Lipschitz condition**

Assume there holds

$$
\left< f(t,y)-f(t,z), y-z \right> \leq \alpha \| y - z \|_2^2
$$

with some $\alpha \in {\mathbb R}$. Then the stability estimate can be improved to

$$
\| y(t) - z(t) \|_2 \leq e^{\alpha (t-t_0)} \| y_0 - z_0 \|_2
$$

This is an essential improvement of the first example above when $a < 0$. The second example satisfies the one-sided Lipschitz condition with $\alpha = 0$.



## Autonomous ODEs

An ODE is called autonomous iff the right hand side does not depend explicitely on $t$, i.e.
$y^\prime(t) = f(y(t))$. Then the time axis can be shifted arbitrarely.

An non-autonomous ODE in ${\mathbb R}^n$ can be transformed to an autonomous one in
${\mathbb R}^{n+1}$ by letting $y_{n+1}(t) := t$ be another unknown state variable, which satisfies
the ODE $y_{n+t}^\prime = 1$ and i.c. $y_{n+1}(t_0) = t_0$.

## An equivalent integral equation

A smooth solution of the ODE satisfies the follong integral (IE) equation, and vice versa:

$$
y(t) = y_0 + \int_{t_0}^t f(s, y(s)) ds \qquad \forall \, t \in [t_0, T]
$$

We call the ODE and the IE formally equivalent. The IE may have a solution which is not differentiable,
and thus not a solution of the ODE (e.g. if $f$ jumps in $\overline t$, then $y$ has a kink, and is not differentiable.

