# Partial differential equations

This chapter gives an introduction to solving partial differential equations (PDEs) using the open source software package [Netgen/NGSolve](https://www.ngsolve.org).

## Installation:

Install a recent Python. Then it should be easy to install NGSolve using

    pip install jupyter numpy scipy matplotlib
    pip install ngsolve
    pip install webgui_jupyter_widgets


To check the installation of NGSolve run in the console:

    python3 -c "import ngsolve; print(ngsolve.__version__)"

Then, open jupyter-notebook, create a new notebook, create and execute a cell with

    from ngsolve import *
    from ngsolve.webgui import Draw
    Draw (unit_cube.shape);


Known issues are
- Use pip3 instead of pip if there is no pip
- If you get an error like `externally-managed-environment`, then either use
virtual environments, or add the flag `--break-system-packages` to the pip command, see [explanation](https://veronneau.org/python-311-pip-and-breaking-system-packages.html)


If local installation does not work, there are alternatives:

- login to a jupyter server from your browser:

  [jupyterhub.cerbsim.com](https://jupyterhub.cerbsim.com) <br>
  user: **ngshub_xx** <br>
  pwd:  **solve!xx** <br>
  with xx number from 01 to 31

  

- run NGSolve online within jupyter-lite:

  [![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://ngsolve.github.io/jupyterlite_ngsolve/lab?path=poisson.ipynb)


  [https://ngsolve.github.io/jupyterlite_ngsolve/lab?path=poisson.ipynb](https://ngsolve.github.io/jupyterlite_ngsolve/lab?path=poisson.ipynb)

  The first time it might take a few minutes to start, and then again to import ngsolve.
  


