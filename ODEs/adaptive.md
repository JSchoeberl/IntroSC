# Automatic step-size control

Some ODEs have local features, which require locally small step sizes. For other regions one can use larger step sizes, 
and save computations. A simple adaptive algorithm is as follows.

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

Compute and plot the trajectory $x(t)$ (recommendation: RK4 method with automatic step-size control). The desired accuracy of $x(T)$ is 1 km.




### Chemical reaction equations: 

molecules of type A and type B react to form molecules of type C and type D.

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

An interesting example is the Belousov-Zhabotinsky reaction (the so called Oregonator) in simplified form is
\begin{eqnarray*}
Br O_3 + Br & \stackrel{k_1}{\longrightarrow } & H Br O_2 + HOBr \\
HBrO_2 + Br & \stackrel{k_2}{\longrightarrow } & 2 HOBr \\
BrO_3 + H Br O_2 & \stackrel{k_3}{\longrightarrow } & 2 HBrO_2 + H_2O \\
HBrO_2 + HBrO_2  & \stackrel{k_4}{\longrightarrow } & HOBr + BrO_3\\
H_2O   & \stackrel{k_5}{\longrightarrow } & Br
\end{eqnarray*}

The concentrations $(c_1, c_2, c_3, c_4, c_5)$ give the densities of ($BrO_3$, $Br$, $HBrO_2$, $HOBr$, $H_2O$), hydrogen is available as much as needed.

This results in the ODE system

\begin{eqnarray*}
c_1^\prime & = & -k_1 c_1 c_2 - k_3 c_1 c_3 + k_4 c_3^2\\
c_2^\prime & = & -k_1 c_1 c_2 - k_2 c_2 c_3 + k_5 c_5 \\
c_3^\prime & = & k_1 c_1 c_2 - k_2 c_2 c_3 + k_3 c_1 c_3 - 2 k_4 c_3^2 \\
c_4^\prime & = & k_1 c_1 c_2 + k_2 c_2 c_3 + k_4 c_3^2 \\
c_5^\prime & = & k_3 c_1 c_3 - k_5 c_5
\end{eqnarray*}

With reaction constants 

$$
k_1 = 1.34, \quad k_2 = 1.6\cdot 10^9, \quad k_3 = 8.0 \cdot 10^3, \quad k_4 = 4.0 \cdot 10^7, \quad k_5 = 1.0, 
$$

Initial conditions are

$$
c_1(0) = 0.05, \quad c_2(0) = 1 \cdot 10^{-4}, \quad c_3(0) = 1 \cdot 10^{-10}, \quad c_4(0) = 0.1, c_5(0) = 1 \cdot 10^{-4}.
$$

Solve the equation with automatic step-size control.
