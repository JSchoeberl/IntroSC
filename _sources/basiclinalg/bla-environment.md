# Setup your development environment


We will develop software in C++ on our own computer, so we have to install developement tools first.


1. [Visual Studio Code](https://code.visualstudio.com) is a popular extensible code editor available for Windows, Linux and macOS.
Please download and install it.


2. Next, we need a compiler. Consult [C/C++ for VS Code](https://code.visualstudio.com/docs/languages/cpp) to install the proper compiler for your system.


3. [CMake](https://cmake.org) is a cross-platform software development tool for building applications, install it. Then, install the *CMake Tools* extension in VS Code.


4. Install *GitLens* VS Code extension module



## Clone a project from github

[Git](https://git-scm.com) is a popular version control system. It is free and open source.

[GitHub](https://github.com) is a platform and cloud-based service for software development and version control using Git.


Look for the github repository
[ASC-bla](https://github.com/TUWien-ASC/ASC-bla)
and git-clone it.

```{image} gitclone.png
:width: 75%
:align: center
```

If everything is setup right, VS Code will first run CMake and then your compiler to build the project.



If you prefer to use command line tools you enter 

    git clone https://github.com/TUWien-ASC/ASC-bla.git


Now we configure and build the project using cmake.
You find a file 'CMakeLists.txt' describing the build structure for the project.

    cd ASC-bla
    mkdir build
    cd build
    cmake ..
    make
    

You should have produced an executable 'test_vector', try it out.



## Extend ASC-bla in team-work

You will now add more functionality to ASC-bla. You want to do it in team-work.
Form teams by 2-4 students. Go back to github, and register. One of the team will *fork* ASC-bla (name it as you like),
and invites colleagues to the project.

Repeat the steps from before: Clone, configure, build, run

Exercise:

 * Extend the library by a `Matrix` class and provide the following operations:
   - Matrix-vector product
   - Matrix-matrix product
   - Transpose matrix function   
   - Some tests for your classes

Allow to choose between row-major and column-major storage:

    enum ORDERING { ColMajor, RowMajor };
    template <typename T, ORDERING ORD>
    class Matrix {
    ...
    }
    

push your changes to your github project using either vs-code functionality, or command line instructions as

    git add demos/demo_matrix.cpp src/matrix.hpp
    git commit -m "added matrix class"
    git push
 
The first command `add` *stages* your files you want to commit to your local version of the repository. With `commit` your local changes are taken over. Finally, with `push` you upload your local commits to the original git repo.


The other team-member can 

    git fetch
    git merge 

the contributions. Experiment with changing the same code regions.


 * Add more advanced features:
   - Inverse matrix using the Gauss-Jordan algorithm: <p>
     Form matrix $M = (A \; I)$, and perform row manipulations (like adding rows, swapping rows and scaling rows) to obtain $\widetilde M = (I \; A^{-1})$.
   - Allow to add `Vector<double>` to a `Vector<std::complex<double>>` <p>
     You can obtain the result-type via
     
         typedef decltype(std::declval<TA>()+std::declval<TB>()) TRES;



