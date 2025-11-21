## Creating documentation

We are using [Jupyter Book](https://jupyterbook.org/en/stable/intro.html) for creating user documentation.
It requires some configuration files, the contents is provided by markdown (.md) and jupyter-notebooks (.ipynb) files.


Install the python packages jupyter-book, and one more to work with github pages:

    pip install -U "jupyter-book<2.0" ghp-import


**Note:** This workflow referes to Jupyter Book 1. The recently (Nov '25) released Jupyter Book 2 requires a different one.

The ASC-bla project includes a little start for the documentation:


Quick steps to setup the documentation are as follows. You don't have to do it for the bla project, files are already there.

    cd docs
    jupyter-book build .


This creates an html-documentation in `_build/html`, which you can now open in the browser.

Adapt the '_config.yml' and '_toc.yml' files for your project.


To publish the documentation within your own bla-project on github run

    ghp-import -n -p -f _build/html




