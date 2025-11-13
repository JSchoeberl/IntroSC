# Overview

We learn how to unleash the full power of our computers. This includes

- **Vectorization**: single instructions can operate on multiple data, for example perform four double precision multiplications simultaneously.
- **Pipelining**: several vector-instructions can be started per clock
  cycle, but it takes several cycles to complete. Start new work before
  results of old work are ready.
- **Caches**: getting data from main memory is slow. Caches keep recently used data in fast but small cache memory. 
- **Parallelization**: most CPUs have several physical cores. Write multi-threaded code to
split the work across multiple cores.


Fork and clone github repo [ASC-HPC](https://github.com/TUWien-ASC/ASC-HPC.git).

This repository contains another github project as a [git submodule](https://www.atlassian.com/git/tutorials/git-submodule),
which you have to initialize using

    git submodule update --init


