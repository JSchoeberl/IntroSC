## Pipelining

Look at the documentation of the `_mm256_fmadd_pd` (fused-multiply-add) intrinsic:
[docu](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=256_fmadd_pd&avxnewtechs=FMA&ig_expand=3101,3101).
It tells us that the instruction has a latency of 4 clock-cycles, but there is a throughput (CPI=clocks per instruction)
of 0.5. This means that we can start two of these instruction in every clock cycle, but we have to wait 4 clock cycles after
a result is available. 

A typical algorithm suffering from latency is computing an inner product:

```cpp
double sum = 0;
for (size_t i = 0; i < n; i++)
   sum += x[i] * y[i];
```

We have to wait until a result is available before we can start the next addition. Thus, one loop takes (at least) 4 clock
cycles. This is called a long dependency chain. We can utilized only $\frac{CPI}{latency} = \frac{1}{8}$ of the possible arithmetic performance.

A vector-vector addition like
```cpp
for (size_t i = 0; i < n; i++)
   x[i] += alpha * y[i];
```
does not suffer from the dependeny chain. The next operation does not depend on the result of the previous one.


A possible optimization for overcoming the dependency chain is unrolling and reordering summation:
```cpp
double sum0 = 0;
double sum1 = 0;
for (size_t i = 0; i+2 <= n; i+=2) {
   sum0 += x[i] * y[i];
   sum1 += x[i+1] * y[i+1];
}
// take care of the leftover
sum = sum0+sum1;
```
We are performing two independent fma - operations per loop, but have to pay the latency only once.
With eight summation variables the full latency bottleneck could be overcome.


However, not only arithmetic operations have a cost, also getting data in and out of the CPU counts.
Modern Intel CPUs can perform per clock cycle at most 4 operations out of

* two simd256 arithmetic operations 
* two simd256 loads and one simd256 store
* four integer operations
* one branch operations

A great source of information are the [manuals by Agner Fog](https://www.agner.org/optimize/#manuals).
Have a look into manual 3. The microarchitecture of Intel, AMD and VIA CPUs: An optimization guide for assembly programmers and compiler makers, Chapter 11: Intel Skylake pipeline.



Since we need more memory transfer operations than arithmetic operations, we cannot gain peak performance
for inner products, and vector addition. E.g. for inner products the ration of memory transfer to arithmetic
operations is $2:1$. This cannot be overcome by reordering summation. If we have to compute inner products
between several pairs of vectors the situation changes:

```cpp
double sum00 = 0;
double sum01 = 0;
double sum10 = 0;
double sum11 = 0;
for (size_t i = 0; i < n; i++) {
   sum00 += x0[i] * y0[i];
   sum01 += x0[i] * y1[i];
   sum10 += x1[i] * y0[i];
   sum11 += x1[i] * y1[i];
}
```

In every loop we load 2 x-values, and 2 y-values, and perform 4 fma operations. Now memory transfer to arithmetic
operations is $1 : 1$. Instead of vectors of scalars, each entry can also be a 256bit `SIMD<double,4>` value.

An important use multiple inner products is a matrix-matrix multiplication. Assuem $A$ is stored row-major, and $B$ is col-major. Then entries $C_{ij}, C_{i+1,j}, C_{i,j+1}, C_{i+1,j+1}$ can be computed by simultaneous inner products of rows $i$ and $i+1$ of $A$, and columns $j$ and $j+1$ of $B$.

The most important case is that $A,B,C$ are stored with the same ordering, say row-major. Then we can multiply two rows of $A$ with $8$ columns (what are two SIMD<double,4>) of $B$ simultaneously:

```cpp
SIMD<double,4> sum00 = 0;
SIMD<double,4> sum01 = 0;
SIMD<double,4> sum10 = 0;
SIMD<double,4> sum11 = 0;
for (size_t i = 0; i < wa; i++) {
   sum00 += rowa0[i] * SIMD<double,4>(colb+i*dist);
   sum01 += rowa0[i] * SIMD<double,4>(colb+i*dist+4);
   sum10 += rowa1[i] * SIMD<double,4>(colb+i*dist);
   sum11 += rowa1[i] * SIMD<double,4>(colb+i*dist+4);
}
```

**Excersicse:**
* Try the examples in the file *demos/simd_timings.cc$ containing such loops.
We measure run-time, and computes GFlop-rate, i.e. billion fma-instructions per second.
* Experiment with different block sizes for simultaneous vector updates. What are your best GFlop-rates ?  
* Look at assembly code, find the inner loops of your functions (enter 'make help' and then 'make demos/simd_timings.s').
* Include the SIMD-classes into your BLA project. Implement a matrix-matrix multiplication using the best loops from simd_timings.cc.


