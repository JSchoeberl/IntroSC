


# The SIMD - class

We do not recommend to write large code using intrinsics, but wrap vector data types and intrinsic functions into high-level C++ classes. This is done in the *simd_xxx.h* header files. For example, the 4 double vector type `__m256d` is wrapped into the template class `SIMD<double,4>`, in file *simd_avx.h*. SIMD is short for **S**ingle **I**nstruction **M**ultiple **D**ata


```cpp
class SIMD<double,4> {
   __m256d val;
  public:
    SIMD (double _val) : val{_mm256_set1_pd(_val)} {};
    SIMD (double v0, double v1, double v2, double v3) : val{_mm256_set_pd(v3,v2,v1,v0)} {  }
    SIMD (__m256d _val) : val{_val} {};    
    SIMD (double const *p) { val = _mm256_loadu_pd(p); }

    auto Val() const { return val; }    
    void Store (double *p) { _mm256_storeu_pd(p, val); }
};

inline SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b)
  { return SIMD<double,4> (_mm256_add_pd(a.Val(), b.Val())); }

```



The Sum4-function is now coded as
```cpp
void Sum4 (double *a, double *b, double *c) {
  SIMD<double,4> va(a);
  SIMD<double,4> vb(b);
  SIMD<double,4> vc = va+vb;
  vc.Store (c);
}
```
Try it on [Compiler explorer](https://compiler-explorer.com/z/9vPaax3fn)

If you are using an Apple silicon processor (manufactured by ARM), your SIMD types are in the file *simd_arm.h*,
and documentation for the intrinsics can be found [here](https://arm-software.github.io/acle/neon_intrinsics/advsimd.html)



### Usage of the SIMD-class

You find examples in the file *demos/demo_simd.cc*.
The SIMD-template is designed such that you can use it for arbitrary (in some range) size.
If it is a native size for your CPU, it will directly use the vector register. If it is larger, the SIMD-type is split into smaller SIMDs.


Creation of SIMD - vector variables:
* Either from simd-size scalar values
* one value which will be broadcasted to the vector
* a pointer from where simd-size values are loaded

```cpp
SIMD<double,4> a(1.,2.,3.,4.);
SIMD<double,4> b(1.0);
double mem[4] = { 10, 11, 12, 13 };
SIMD<double,4> c(&mem[0]);
```

Output:
```cpp
std::cout << "a = " << a << std::endl;
```

You can perform arithmetic operations with this vector variables (not everything is implemented,
but can be added easily following the same pattern):
```cpp
auto d = a*b+c;
```
The operations are applied to each scalar component of the vector.

Then you can store the result to memory:
```cpp
double res[4];
d.Store(&res[0]);
```

### Masking operations

We can generate boolean-SIMDs like this:
```cpp
SIMD<int64_t,4> sequ = IndexSequence<int64_t, 4>();
cout << "sequ = " << sequ << endl;
SIMD<mask64,4> mask = (2 >= sequ);
cout << "2 >= " << sequ << " = " << mask << endl;
```
The output is as follows:
```cpp
sequ = 0, 1, 2, 3
2 >= 0, 1, 2, 3 = t, t, t, f
```

The mask can be used to select component-wise between two simd vectors:
```cpp
auto res = Select(mask, a, b);
```
For the mask `(t,t,t,f)`, the result `res` is `(a0,a1,a2,b3)`.

A mask can also be used to load or store only a subset of vector components.
```cpp
c.Store(&mem[0], mask);
```

Such a mask is useful to treat the left-over of a long vector operations. Let's copy
vector x to vector y. The bulk of the vector is processed using SIMD<4> instructions.
For the $n\%4$ elements at the end one could start a scalar loop, or use masking:
```cpp
size_t i = 0;
for ( ; i+4 <= n; i += 4)
   SIMD<double,4> (px+i).Store(py+i);
auto mask = (n%4) > IndexSequence<int64_t, 4>();
SIMD<double,4>(px+i, mask).Store(py+i, mask);
```

### Horizontal operations, shuffling

SIMD - operations work very well when data is aligned properly. Sometimes one has to move data
around the individual components, what is less fun (and less efficient). An example is summing
all entries of a
simd-vector, called horizontal sum. This functions forms the horizontal sums of two given vectors:
```cpp
SIMD<double,2> HSum (SIMD<double,4> sd1, SIMD<double,4> sd2) {
  __m256d hv = _mm256_hadd_pd(sd1.Val(), sd2.Val());
  __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
  return SIMD<double,2>(_mm_cvtsd_f64 (hv2), _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
```
Lookup the documentation of the used intrinsics.

**Exercise:**

* Implement missing arithmetic operations for the SIMD template class

* (Intel CPUs only) Implement `SIMD<T,2>` classes similar to the 4-wide SIMDs. Use intrinsics to construct/split the `SIMD<T,4>` from/to `SIMD<T,2>` types.

* (ARM CPUs only) Implement missing comparison operators

* Implement a $4\times4$ - Transpose operations: The 4 input SIMD-vectors `ai` contain the rows of a $4\times4$-matrix. The output vectors `bi` should have the rows of the transpose matrix:
```
void Transpose (SIMD<double,4> a0, SIMD<double,4> a1, SIMD<double,4> a2, SIMD<double,4> a3,
                SIMD<double,4> &b0, SIMD<double,4> &b1, SIMD<double,4> &b2, SIMD<double,4> &b3);
```
Useful functions are `_mm256_unpacklo_pd`, `_mm256_unpackhi_pd`,
`_mm256_insertf128_pd`, and `_mm256_extractf128_pd`.





