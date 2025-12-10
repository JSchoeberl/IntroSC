# Automatic step-size control

Some ODEs have local features, which require locally small step sizes. For other regions on can use larger step sizes, 
and save computations. A simple adaptive algorithm is as follows.

To estimate the local error, we apply two different time stepping methods, and take the difference as error estimation. 
The error in one time-step should be about $\tau \times \text{tol}$.


```
while t < Tend
    y1 = y;
    method1.doStep(tau, y1);
    y2 = y;
    method2.doStep(tau, y2);
    double errest = norm(y2-y1);

    if errest < tau*tol
        t = t + tau;
        y = y1;  // accept step
        tau = tau * 1.5;  // increase step size
    else
        tau = tau * 0.5;  // decrease step size       
```

A better way to increase the time-step is

$$
  \tau = \sqrt[p+1]{\frac{0.5 * tol}{errest}} \; \tau
$$

If Newton's method did not converge, we also try with a smaller time-step:
```code
try {
  method1.doStep(tau,y1);
  method2.doStep(tau,y2);
}
catch (std::exception & e) {. // Nowton did not converge
  tau *= 0.5;
  continue;
}
```




## Examples

### Arenstorf-Orbit 

*vgl. Aufgabe 4.4 aus Numerische Mathematik II, P. Deuflhard, F. Bornemann*
Die Differentialgleichungen der Bewegung eines Satelliten um das System Erde-Mond lautet in den Koordinaten $x =(x_1,x_2)$ des mitrotierenden Schwerpunktsystems
\begin{align*} 
x_1'' & = x_1 + 2 x_2' - \hat{\mu} \frac{x_1 + \mu}{N_1}- \mu \frac{x_1-\hat{\mu}}{N_2} \\
x_2'' &= x_2 - 2 x_1' - \hat{\mu} \frac{x_2}{N_1} - \mu \frac{x_2}{N_2}
\end{align*}
mit  den Abkürzungen
\begin{equation*}
N_1 = ( (x_1 + \mu)^2 + x_2^2)^{3/2}, \quad N_2 = ( (x_1 - \hat{\mu})^2 + x_2^2 )^{3/2}
\end{equation*}
und den Daten
\begin{equation*}
\mu = 0.012277471 , \quad \hat{\mu} = 1  - \mu.
\end{equation*}
Dabei ist $\mu$ das Verhältnis der Mondmasse zur Masse der Gesamtsystems. L\"angeneinheit f\"ur die euklidische Norm $|x|$ ist die mittlere Entfernung Erde-Mond (ca. $384000 \, km$), 
Zeiteinheit ein Monat. 
Die Anfangswerte
\begin{equation*}
x_1(0) = 0.994, \quad x_1'(0) = 0, \quad x_2(0) = 0, \quad x_2'(0) = -2.0015851063790825
\end{equation*}
sind so gewählt, dass sich der kleeblattförmige [Arenstorf-Orbit ergibt](http://en.wikipedia.org/wiki/Richard_Arenstorf). Die Periode dieses Orbits beträgt
\begin{equation*}
T = 17.0652166.
\end{equation*}

Berechnen und plotten Sie die Bahnkurve $x(t)$ (Empfehlung: RK4 Verfahren und automatische Schrittweitensteuerung). Es wird eine Genauigkeit von $x(T)$ von 1 km gewünscht.




### Chemical reaction equations: 

molecules of type A and type B react to form molecules of type C and type D.

$$
A + B \stackrel{k}{\longrightarrow } C + D
$$

The probability of the reaction is $k c_A c_B$, where $c_A, c_B, ...$ are concentrations of molecules.
This leads to differential equations:

\begin{eqnarray*}
c_A^\prime & = & - k c_A c_B \\
c_B^\prime & = & - k c_A c_B \\
c_C^\prime & = &   k c_A c_B \\ 
c_D^\prime & = &  k c_A c_B. 
\end{eqnarray*}

An interesting example is the  Zhabotinski-Belousov-Reaction (the so called Oregonator) in simplified form is
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
