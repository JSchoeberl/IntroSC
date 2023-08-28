# Basic Linear Algebra

We are developing a library for basic linear algebra operations. For this we
* are implementing classes for vectors and matrices in modern C++
* make them available to Python for simple usability

We are practicing software development in teams using git.


## Clone a prepared project from github

[Git](https://git-scm.com) is a popular version control system. It is free and open source.

[GitHub](https://github.com) is a platform and cloud-based service for software development and version control using Git.


Prerequisites:
* install a C++ environment on your computer (e.g. Visual Studio Code)
* install git (available within vs-code as ...)
* install cmake (available within vs-code as ...)

look for project draft 
[ASC-bla](https://github.com/JSchoeberl/ASC-bla)

and git-clone it. Either
* use vscode ...
* git clone <copy git link>
* on Windows: ...


Now we configure and build the project using cmake.
You find a file 'CMakeLists.txt' describing the build structure for the project.

Again, *vscode* is doing the work for you. If you use a terminal, go to the ASC-bla directory and run:

    mkdir build
    cd build
    cmake ..
    make
    

You should now be able to run the executeable.

- [Refresh your C++](bla-explainCpp.md)


## Extend ASC-bla in team-work

You will now add more functionality to ASC-bla. You want to do it in team-work.
Form teams by 2-4 students. Go back to github. One of the team will *fork* ASC-bla (name it as you like),
and invites collaegues to the project.

Repeat the steps from before: Clone, configure, build, run

Homework:

 * Extend the library by a **Matrix** class, and provide the following operations:
   - Matrix-Vector product
   - Transpose matrix function
   - Some tests for your classes

Commit and push your changes to your github project using either vs-code functionality, or command line instructions as

    git add matrix.cc matrix.h
    git commit -m "added matrix class"
    git push
 

The other team-member can 

    git fetch
    git merge 

the contributions. Experiment with changing the same code reagions. 

## Preparing documentation

The Together with 
Read the Docs

Together with the 







