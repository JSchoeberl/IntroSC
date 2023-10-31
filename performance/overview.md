# Overview

We learn how to unleash the full power of our computers. This include

- **Vectorization**: single instructions can operate with multiple data, for example perform four double precision multiplications simulteneously
- **Pipelining**: several vector-instructions can be started per clock cycle, but it takes several cycles to complete. Start new work befor old work is ready
- **Caches**: getting data from main memory is expensive. Caches keep recently used data in fast but small cache memory. 
- **Parallelization**: most CPUs have several physical cores. Write multi-threaded code to split work to multiple cores


Fork and clone github repo [ASC-HPC](https://github.com/TUWien-ASC/ASC-HPC.git).

This repository contains another github project as a submodule,
which you have to initialize using

    git submodule update --init


