## Creating documentation

We use [ReadTheDocs](https://docs.readthedocs.io) based on [Sphinx](https://www.sphinx-doc.org) for creating documentation.

Install the python packages *sphinx* and the rtd-theme using

    pip install sphinx sphinx_rtd_theme

Quick steps to setup the documentation are as follows. You don't have to do it for the bla project, files are already there.

    mkdir docs
    cd docs
    sphinx-quickstart

Now, and after changing contents in the file `index.rst` run 

    make html


In the configuration file `conf.py` you can select the theme, for example:

    html_theme = 'sphinx_rtd_theme'


To publish your documentation online register at [readthedocs.org](https://readthedocs.org)
with your github account.
On your dash-board you can see your (public!) github projects.

Import your bla-project, and give a unique name for the documentation.

Put the 'docs passing' batch into your github-readme file.





