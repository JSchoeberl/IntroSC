## Clone a project from github

[Git](https://git-scm.com) is a popular version control system. It is free and open source.

[GitHub](https://github.com) is a platform and cloud-based service for software development and version control using Git.


Prerequisites:
* install a C++ environment on your computer (e.g. [Visual Studio Code](https://code.visualstudio.com))
* install git (available within vs-code as *Git Extension Pack*)
* install cmake (available within vs-code as *CMake Tools*)

look for project draft 
[ASC-bla](https://github.com/TUWien-ASC/ASC-bla)

and git-clone it. Either
* use vscode 'source control -> Clone Repository'
* git clone \<copy git link\>
* on Windows: ...


Now we configure and build the project using cmake.
You find a file 'CMakeLists.txt' describing the build structure for the project.

Again, *vscode* is doing the work for you. If you use a terminal, go to the ASC-bla directory and run:

    mkdir build
    cd build
    cmake ..
    make
    

You should have produced an executable 'test_vector', try it out.



### Extend ASC-bla in team-work

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

    git add tests/test_matrix.cc src/matrix.h
    git commit -m "added matrix class"
    git push
 
The first command `add` *stages* your files you want to commit to your local version of the repository. With `commit` your local changes are taken over. Finally, with `push` you upload your local commits to the original git repo.


The other team-member can 

    git fetch
    git merge 

the contributions. Experiment with changing the same code regions.


 * Add more advanced features:
   - Inverse matrix: <p>
     Form matrix $M = (A \; I)$, and perform row manipulations to obtain $\widetilde M = (I \; A^{-1})$.
   - Allow to add `Vector<double>` to a `Vector<std::complex<double>>`

         typedef decltype(std::declval<TA>()+std::declval<TB>()) TRES;



