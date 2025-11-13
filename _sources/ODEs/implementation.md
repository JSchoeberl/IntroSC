# Implementation

We give an implementation of time-stepping methods.

* We introduce the class `NonlinearFunction`. It represents the right hand side of the ODE. It allows to compose functios to build up complex expressions from simle building blocks.

* We introduce a class hierarchy of `TimerStepper`s for implementing different methods
