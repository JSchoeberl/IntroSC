## Clone a project from github

[Git](https://git-scm.com) is a popular version control system. It is free and open source.

[GitHub](https://github.com) is a platform and cloud-based service for software development and version control using Git.


Prerequisites:
* install a C++ environment on your computer (e.g. Visual Studio Code)
* install git (available within vs-code as ...)
* install cmake (available within vs-code as ...)

look for project draft 
[ASC-bla](https://github.com/TUWien-ASC/ASC-bla)

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
    

You should now be able to run the executable.



### Extend ASC-bla in team-work

You will now add more functionality to ASC-bla. You want to do it in team-work.
Form teams by 2-4 students. Go back to github. One of the team will *fork* ASC-bla (name it as you like),
and invites colleagues to the project.

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

the contributions. Experiment with changing the same code regions. 