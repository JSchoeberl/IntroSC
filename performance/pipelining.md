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



