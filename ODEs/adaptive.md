# Automatic step-size control

Some ODEs have local features, which require locally small step sizes. For other regions one can use larger step sizes,
and save computations. A simple adaptive algorithm is as follows.

## Estimating the local error

To estimate the local error, we apply two different time stepping methods, and take the difference as error estimation.
Method 1 has order $p$, method 2 is more accurate. Then the difference estimates the local error of method 1, which behaves like $C \tau^{p+1}$.
The error in one time-step should be about $\tau \times \text{tol}$.


```
while t < Tend
    tau = min(tau, Tend-t);
    y1 = y;
    method1.doStep(tau, y1);
    y2 = y;
    method2.doStep(tau, y2);
    double errest = norm(y2-y1);

    if errest < tau*tol
        t = t + tau;
        y = y2;  // accept step, continue with the more accurate solution
        tau = tau * 1.5;  // increase step size
    else
        tau = tau * 0.5;  // decrease step size
```

A better way to increase the time-step is

$$
  \tau_{\text{new}} = \sqrt[p]{\frac{0.5 \, \tau \, \text{tol}}{\text{errest}}} \; \tau,
$$

which aims at $\text{errest} \approx 0.5 \, \tau_{\text{new}} \, \text{tol}$ in the next step.

If Newton's method did not converge, we also try with a smaller time-step:
```cpp
try {
  method1.doStep(tau,y1);
  method2.doStep(tau,y2);
}
catch (std::exception & e) { // Newton did not converge
  tau *= 0.5;
  continue;
}
```

### In practice

Production codes use a few more ingredients:

* The tolerance is often posed per step, i.e. accept if $\text{errest} \leq \text{tol}$. Then
  $\tau_{\text{new}} = \tau \, \big( \text{fac} \cdot \text{tol} / \text{errest} \big)^{1/(p+1)}$, with a safety factor like $\text{fac} = 0.9$.

* The change of the time-step is limited, for example $0.2 \, \tau \leq \tau_{\text{new}} \leq 5 \, \tau$.

* If the components of $y$ have very different magnitudes (in the chemical example below from $10^{-11}$ to $10^{-1}$),
  the Euclidean norm only sees the large components. One uses absolute and relative tolerances, and the scaled norm

  $$
  \text{errest} = \Big( \frac{1}{d} \sum_{k=1}^d \Big( \frac{y_{2,k} - y_{1,k}}{\text{atol} + \text{rtol} \, |y_k|} \Big)^2 \Big)^{1/2},
  $$

  which is compared with $\text{tol} = 1$.


## Embedded Runge-Kutta methods

Computing two independent methods is expensive. Embedded Runge-Kutta methods share all stages $k^j$,
i.e. the $a$ and $c$ coefficients, and use two different weight vectors $b$ and $\hat b$:

$$
y_{n+1} = y_n + \tau \sum_j b_j k^j, \qquad \hat y_{n+1} = y_n + \tau \sum_j \hat b_j k^j.
$$

The weights are chosen such that $y_{n+1}$ has order $p+1$, and $\hat y_{n+1}$ has order $p$.
The error estimate

$$
\text{errest} = \| y_{n+1} - \hat y_{n+1} \| = \tau \, \Big\| \sum_j (b_j - \hat b_j) k^j \Big\|
$$

comes for free, no additional evaluations of $f$ are needed. The estimate is for the lower order
solution $\hat y_{n+1}$, so $p$ is the exponent in the step-size formula, but we continue with the more accurate $y_{n+1}$.
The Butcher tableau is extended by one row for $\hat b$.

The simplest example is the Heun-Euler pair of orders 2 and 1:

$$
\begin{array}{c|cc}
0 & & \\
1 & 1 & \\
\hline
 & \frac{1}{2} & \frac{1}{2} \\
 & 1 & 0
\end{array}
$$

The Bogacki-Shampine pair has orders 3 and 2:

$$
\begin{array}{c|cccc}
0 & & & & \\
\frac{1}{2} & \frac{1}{2} & & & \\
\frac{3}{4} & 0 & \frac{3}{4} & & \\
1 & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} & \\
\hline
 & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} & 0 \\
 & \frac{7}{24} & \frac{1}{4} & \frac{1}{3} & \frac{1}{8}
\end{array}
$$

The last row of $a$ coincides with $b$. Thus the last stage is $k^3 = f(y_{n+1})$, which is also the first stage
of the next time-step. This trick is called *first same as last* (FSAL), and the method needs only three evaluations of $f$ per step.

The most popular embedded method is the 7-stage Dormand-Prince pair of orders 5 and 4, also using FSAL.
It is the default method `ode45` in Matlab, and `RK45` in scipy.


## Stiff problems

Adaptive explicit methods do not fail for stiff problems, but they become very expensive:
When $\tau \lambda$ leaves the domain of stability, the numerical solution starts to grow, the error estimator
detects this, and reduces the time-step. Thus the time-step is determined by stability, and not by accuracy,
even when the solution is smooth.

For the Oregonator example below, an explicit adaptive 3(2) method needs about one million steps on the interval $[0, 50]$,
while an implicit Radau IIA method needs about 300 steps (both computed with scipy, $\text{rtol} = 10^{-6}$).

For stiff problems we use a pair of implicit methods, preferably both $L$-stable. An example is the implicit
Euler method ($p = 1$) for the error estimate, together with the two-stage Radau IIA method (order 3).


## Examples

### Arenstorf-Orbit

*cf. Exercise 4.4 in P. Deuflhard, F. Bornemann, Numerische Mathematik II*

The equations of motion of a satellite in the earth-moon system, in coordinates $x =(x_1,x_2)$ of the co-rotating barycentric frame, are
\begin{align*}
x_1'' & = x_1 + 2 x_2' - \hat{\mu} \frac{x_1 + \mu}{N_1}- \mu \frac{x_1-\hat{\mu}}{N_2} \\
x_2'' &= x_2 - 2 x_1' - \hat{\mu} \frac{x_2}{N_1} - \mu \frac{x_2}{N_2}
\end{align*}
with the abbreviations
\begin{equation*}
N_1 = ( (x_1 + \mu)^2 + x_2^2)^{3/2}, \quad N_2 = ( (x_1 - \hat{\mu})^2 + x_2^2 )^{3/2}
\end{equation*}
and the data
\begin{equation*}
\mu = 0.012277471 , \quad \hat{\mu} = 1  - \mu.
\end{equation*}
Here, $\mu$ is the ratio of the moon's mass to the mass of the total system. The unit of length for the Euclidean norm $|x|$ is the mean earth-moon distance (about $384000 \, km$),
the unit of time is one month.
The initial values
\begin{equation*}
x_1(0) = 0.994, \quad x_1'(0) = 0, \quad x_2(0) = 0, \quad x_2'(0) = -2.0015851063790825
\end{equation*}
are chosen such that the cloverleaf-shaped [Arenstorf orbit](http://en.wikipedia.org/wiki/Richard_Arenstorf) is obtained. The period of this orbit is
\begin{equation*}
T = 17.0652166.
\end{equation*}

The solution is non-stiff, but it has phases of very fast motion close to the earth and the moon.


### Chemical reaction equations

Molecules of type A and type B react to form molecules of type C and type D:

$$
A + B \stackrel{k}{\longrightarrow } C + D
$$

The reaction rate is $k c_A c_B$, where $c_A, c_B, ...$ are concentrations of molecules.
This leads to differential equations:

\begin{eqnarray*}
c_A^\prime & = & - k c_A c_B \\
c_B^\prime & = & - k c_A c_B \\
c_C^\prime & = &   k c_A c_B \\
c_D^\prime & = &  k c_A c_B.
\end{eqnarray*}

An interesting example is the Belousov-Zhabotinsky reaction. Its simplified model, the so called Oregonator, is

\begin{eqnarray*}
\mathrm{BrO_3^- + Br^-} & \stackrel{k_1}{\longrightarrow } & \mathrm{HBrO_2 + HOBr} \\
\mathrm{HBrO_2 + Br^-} & \stackrel{k_2}{\longrightarrow } & \mathrm{2 \, HOBr} \\
\mathrm{BrO_3^- + HBrO_2} & \stackrel{k_3}{\longrightarrow } & \mathrm{2 \, HBrO_2 + Ce^{4+}} \\
\mathrm{2 \, HBrO_2} & \stackrel{k_4}{\longrightarrow } & \mathrm{HOBr + BrO_3^-} \\
\mathrm{Ce^{4+}} & \stackrel{k_5}{\longrightarrow } & \mathrm{Br^-}
\end{eqnarray*}

The concentrations $(c_1, c_2, c_3, c_4, c_5)$ give the densities of ($\mathrm{BrO_3^-}$, $\mathrm{Br^-}$, $\mathrm{HBrO_2}$, $\mathrm{HOBr}$, $\mathrm{Ce^{4+}}$).
The cerium ions act as catalyst, the last reaction summarizes the oxidation of an organic substance (malonic acid), which
releases bromide. Hydrogen ions and the organic substance are available as much as needed.

This results in the ODE system

\begin{eqnarray*}
c_1^\prime & = & -k_1 c_1 c_2 - k_3 c_1 c_3 + k_4 c_3^2\\
c_2^\prime & = & -k_1 c_1 c_2 - k_2 c_2 c_3 + k_5 c_5 \\
c_3^\prime & = & k_1 c_1 c_2 - k_2 c_2 c_3 + k_3 c_1 c_3 - 2 k_4 c_3^2 \\
c_4^\prime & = & k_1 c_1 c_2 + 2 k_2 c_2 c_3 + k_4 c_3^2 \\
c_5^\prime & = & k_3 c_1 c_3 - k_5 c_5
\end{eqnarray*}

With reaction constants

$$
k_1 = 1.34, \quad k_2 = 1.6\cdot 10^9, \quad k_3 = 8.0 \cdot 10^3, \quad k_4 = 4.0 \cdot 10^7, \quad k_5 = 1.0,
$$

Initial conditions are

$$
c_1(0) = 0.05, \quad c_2(0) = 1 \cdot 10^{-4}, \quad c_3(0) = 1 \cdot 10^{-10}, \quad c_4(0) = 0.1, \quad c_5(0) = 1 \cdot 10^{-4}.
$$

The concentrations oscillate with a period of roughly 60 time units. The eigenvalues of the Jacobi-matrix reach $-4 \cdot 10^5$, the system is stiff.


## Exercises

* Implement the adaptive time-stepping loop for two `TimeStepper` objects, including the handling of non-converging Newton iterations.

* Extend your `ExplicitRungeKutta` class to an embedded method with weights $b$ and $\hat b$, which also provides the error estimate.
  Implement the Bogacki-Shampine pair, and verify the orders of $b$ and $\hat b$ with constant time-steps.

* Compute and plot the Arenstorf orbit $x(t)$ with automatic step-size control. The desired accuracy of $x(T)$ is 1 km.
  Plot the time-step $\tau$ over time. How many steps does the classical $RK_4$ method with constant time-step need for the same accuracy?

* Solve the Oregonator equations on the interval $[0, 300]$ with automatic step-size control, using the scaled error norm.
  Plot the concentrations in logarithmic scale. Compare the number of time-steps of the explicit Bogacki-Shampine pair
  with the implicit pair of implicit Euler and Radau IIA.
