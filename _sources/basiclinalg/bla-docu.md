## Creating documentation

We use [ReadTheDocs](https://docs.readthedocs.io) based on [Sphinx](https://www.sphinx-doc.org) for creating documentation.

Install the python packages *sphinx* and the rtd-theme using

    pip install sphinx sphinx_rtd_theme

Quick steps to setup the documentation are as follows. You don't have to do it for the bla projects, files are already there. Only run make after changing the documentation files. 

    mkdir docs
    cd docs
    sphinx-quickstart
    make


in `conf.py` you can select the theme, for example:

    html_theme = 'sphinx_rtd_theme'


Register to [readthedocs.org](readthedocs.org) with your github account.
On your dash-board you can see your (public!) github projects. Import your
bla-project. 




