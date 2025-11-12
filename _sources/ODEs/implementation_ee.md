## Implementing the explicit Euler method

We want to implement a time-stepping method for solving the autonomous ODE

$$
y\prime (t) = f(y(t)) \qquad \forall t \in [0, T], \qquad y(0) = y_0
$$

using the explicit Euler method:

$$
y_{i+1} = y_i + \tau f(y_i) \qquad 0 \leq i < n
$$

with $n = T/\tau$.


A simple implementation of an explicit Euler method may look like that;


```cpp
void SolveODE_EE(double tend, int steps,
                 VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                 std::function<void(double,VectorView<double>)> callback = nullptr)
{
  double dt = tend/steps;
  Vector<double> fy(y.size());

  double t = 0;
  for (int i = 0; i < steps; i++)
    {
      rhs -> evaluate (y, fy);
      y += dt * fy;
      t += dt;
      if (callback) callback(t, y);
    }
}
```




## Using the time-stepping method

A mass attached to a spring is described by the ODE

$$
m y^{\prime \prime}(t) = -k y(t)
$$

where $m$ is mass, $k$ is the stiffness of the spring, and $y(t)$ is the
displacement of the mass. The equation comes from Newton's law

> force = mass $\times$ acceleration 

We replace the second-order equation with a first order system. For mathematical simplicity we set $k = m = 1$.
Then we can define the right-hand-side as a `NonlinearFunction`. The derivative is needed for the Newton solver:

```cpp
class MassSpring : public NonlinearFunction
{
  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};
```

Finally, We start the time-stepper with the time interval, number of steps, initial values,
right-hand-side function, and a call-back function called at the end of every time-step:
```cpp
double tend = 4*M_PI;
int steps = 100;
Vector<> y { 1, 0 };
auto rhs = make_shared<MassSpring>();

SolveODE_IE(tend, steps, y, rhs,
            [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
```
This example is provided in `demos/test_ode.cc`.





