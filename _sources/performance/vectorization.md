# Vectorization

This function adds four doubles:
```cpp
double a[4], b[4], c[4];
void Sum4() {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
  c[3] = a[3] + b[3];
}
```

Compiling it for modern intel CPUs gives this assembly code:
[Compiler explorer](https://compiler-explorer.com/z/hjvM1zo9P)


```asm
vmovapd ymm0, ymmword ptr [rip + a]
vaddpd  ymm0, ymm0, ymmword ptr [rip + b]
vmovapd ymmword ptr [rip + c], ymm0
```
The first instruction loads four doubles from variable `a` into the 256-bit register `ymm0`, the second one (`vaddpd`) adds four doubles from variable `b` to register `ymm0`, and the third one stores the result to variable `c`.
The compiler knows what kind of instructions are supported by the CPU.

`ymm` are CPU registers of width 256 bits, they can hold four `double`s (double precision floating point numbers), or eight `float`s (single precision floating point numbers). These 256bit registers are available in CPUs supporting the [AVX extension](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions).
The older `xmm` registers can hold two doubles, and `zmm` registers can even hold 8 doubles, they are available on server CPUs, and only within a few consumer CPUs.

Now, try this similar function (copy into compiler explorer):
```cpp
void Sum4(double *a, double *b, double *c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
  c[3] = a[3] + b[3];
}
```
why can't the compiler generate similar short code ?

Modern compilers can do a lot of auto-vectorization, but often the outcome is more predictable if we do the vectorization
by hand. For that, compiler provide vector data types, and so called **intrinsic** functions, which look like C functions:
```cpp
#include <immintrin.h>
void Sum4(double *a, double *b, double *c) {
   __m256d va = _mm256_loadu_pd(a);
   __m256d vb = _mm256_loadu_pd(b);
   __m256d vc = _mm256_add_pd(va, vb);
   _mm256_storeu_pd(c, vc);
}
```
On Intel (and compatible) CPUs we include the the *immintrin.h* header. This provides the data-type `__mm256d` for a vector of four doubles (which are 256bits). The intrinsic function `_mm256_loadu_pd` loads such a 4-vector starting from the address given by a double pointer. The intrinsic function `_mm256_add_pd` adds two 4-vectors, and finally ` _mm256_storeu_pd` stores the 4-vector in memory, starting at the address given by a double pointer. Try this code in compiler explorer !


Available intrinsics are found in the
[Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html).
Lookup the documentation of the functions used above!

These intrinsics are supported by all major compilers (gcc, clang, msvc, icc, ...), however with slightly different behaviour. For example, gcc has arithmetic operations for type `__m256d` predefined, in msvc the intrinsic function must be called.


Vectorization is like
```{image} skilift.jpg
:width: 70%
:align: center
```