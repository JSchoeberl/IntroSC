# Mechanical Systems

We consider a system of point masses $m_i$ at time-dependent positions $x_i(t)$.

There is potential energy $U(x)$ stored in that system. Examples are gravitational energy

$$
\sum m_i g x_{i,z},
$$

or energy stored in springs $k_{ij}$ connecting masses $m_i$ and $m_j$:

$$
\sum_{k_{ij}} \frac{1}{2} k_{ij} ( |x_i - x_j | - l_{ij} )^2,
$$

where $l_{ij}$ is the length of the spring without load, and $|x_i - x_j|$ is the Euclidean distance.


Taking the negative derivative

$$
F_i := -\frac{\partial U}{\partial x_i}(x)
$$

gives the force onto the mass.

In case of gravity it is $-m_i g e_z$, in case of springs it is

$$
-\sum_j k_{ij} ( |x_i - x_j | - l_{ij} ) \frac{x_i - x_j}{|x_i - x_j|}
$$

i.e. spring constant multiplid by elongation, in the direction of the vector $x_j-x_i$.


By Newton's principle, the force leads to acceleration:

$$
m_i \ddot x_i = -\frac{\partial U}{\partial x_i}
$$

Now we denote the derivative w.r.t. time $t$ via $\dot {x}$.

## The Newmark method

We start with reducing the second order equation to a first order system. Therefore, we introduce velocity $v_i = \dot x_i$:

$$
\begin{array}{rcl}
\dot x_i & = & v_i \\
m_i \dot v_i & = & -F_i(x)
\end{array}
$$

Applying the Crank-Nicolson time-stepping scheme we obtain

$$
\begin{array}{rcl}
x_{n+1} - x_n & = & \frac{h}{2} (v_n + v_{n+1} ) \\
M (v_{n+1} - v_n) & = & \frac{h}{2} (F(x_n) + F(x_{n+1}))
\end{array}
$$

By renaming $a_n := M^{-1} F(x_n)$, we can substitute 

$$
\begin{array}{rcl}
v_{n+1} & = & v_n + \frac{h}{2} (a_n + a_{n+1}) \\
x_{n+1} & = & x_n + \frac{h}{2} (v_n + v_{n+1}) = x_n + h v_n + \frac{h^2}{4} (a_n + a_{n+1} )
\end{array}
$$

and obtain the equation for the new $a_{n+1}$:

$$
M a_{n+1} = F \big( x_n + h v_n + \tfrac{h^2}{4} (a_n + a_{n+1} ) \big)
$$


## Generalized-$\alpha$ method

Although the Newmark method preserves energy exactly for linear ODEs, it leads to
instabilities for non-linear equations. For this reason, the following generalization is
usually used in mechanical simulation codes:

By adding a bunch of parameters, one ends up with the Generalied $\alpha$ method:

$$
M a_{n+1-\alpha_m} = F(x_{n+1-\alpha_f})
$$

with

$$
\begin{array}{rcl}
x_{n+1} & = & x_n + h v_n  + h^2 ( ( \frac{1}{2} - \beta) a_n + \beta a_{n+1} ) \\
v_{n+1} & = & v_n + h ( (1-\gamma) a_n + \gamma a_{n+1} ) \\
x_{n+1-\alpha_f} & = & (1-\alpha_f) x_{n+1} + \alpha_f x_n \\
a_{n+1-\alpha_m} & = & (1-\alpha_m) a_{n+1} + \alpha_m a_n 
\end{array}
$$

Depending of one free parameter $\rho^\infty$, parameters tuned for optimal accuracy
and stabilty are obtained by

$$
\begin{array}{rcl}
\alpha_m & = & \frac{2 \rho^\infty - 1}{\rho^\infty + 1} \\
\alpha_f & = & \frac{\rho^\infty}{\rho^\infty + 1} \\
\beta & = & \frac{1}{4} (1 - \alpha_m + \alpha_f)^2 \\
\gamma & = & \frac{1}{2} - \alpha_m + \alpha_f
\end{array}
$$

The parameter $\rho^\infty$ specifies the damping for high frequency functions.
The reference is J. Chung and G.M. Hulbert *A Time Integration Algorithm for Structural Dynamics With Improved Numerical Dissipation: The Generalized-α Method*, J. Appl. Mech., 1993, [see](https://asmedigitalcollection.asme.org/appliedmechanics/article/60/2/371/423023/A-Time-Integration-Algorithm-for-Structural?casa_token=PkMuU1hl_4AAAAAA:btyl_IyLpSdrI9q5LjNYSDLzOECUIoT1b3IRqcF86A773ZNWRnyuN94Tizg4kQUaVw2gZftWhg)
An easy to read derivation is [here](https://miaodi.github.io/finite%20element%20method/newmark-generalized/).




## Systems with constraints

Think of a pendulum, where a mass $m$ is attached to a (mass-less) stick of length $l$.
One can describe the system by the angle $\alpha$ between the stick and the vertical axis.
By some high-school mathematics, one splits the total gravitational force into the direction of
the stick, and the orthogonal complement, and ends up with the equation for the angle $\alpha$:

$$
\ddot \alpha = -\frac{g}{l} \sin \alpha
$$

However, such a direct modeling with angles (so called generalized coordinates) becomes
difficult for more complex systems. A simpler possibility is to pose the length constraint

$$
0 = g(x) := | x - x_0 | - l,
$$

where $x_0$ is the anchor point of the pendulum, and define the Lagrange function

$$
L(x, \lambda) = -U(x) + \left< \lambda , g(x) \right>,
$$

with $U(x) = m g x_z$, and pose the equation

$$
m_i \ddot x_i & = & \frac{\partial}{\partial x_i} L(x, \lambda) \\
 0 & = & \nabla_\lambda L(x, \lambda)
$$

The physical meaning of the Lagrange parameter $\lambda$ is the longitudinal force in the pendulum. This is an ordinary differential equation, with algebraic constraints,
a so called differential-algebraic equation (DAE).

One can diretly use the generalized-$\alpha$ method for this DAE. The right hand side is the gradient of the Lagrange function, the mass matrix on the left hand side is extended by $0$ for the Lagrange parameters:

$$
\left( \begin{array}{cc}
 M & 0 \\
 0 & 0
 \end{array} \right)
\left( \begin{array}{cc}
   \ddot x \\
   \ddot \lambda
 \end{array} \right) =
 \nabla L (x, \lambda)
$$




