# Expression templates

Think of a code like


```cpp
z = x + 3*y;
```

where `x,y,z` are vectors. What happens ? First, the `operator* (double, Vector)` is called, it generates a temporary vector for 3*y, and then the `operator+ (Vector, Vector)` is called to compute a new vector for the sum. Finally, the later vector is copied into the memory of `z`.

Writing an old-style C code
```cpp
for (size_t i = 0; i < z.size(); i++)
  z(i) = x(i) + 3*y(i);
```
seems to be much more efficient, no temporary objects are created, and the result values are stored at the right position immediately.


**Expression templates** provide a way to combine the elegance of the first version, with the performance of the explicit loop.

The trick is that `x+3*y` does not return vectors, but it returns the logic information of the expression, like *I am the sum of a vector and a second vector scaled by the scalar 3*. Such expressions can be built by the compiler using templates, something like the type

```cpp
SumExpr< Vector , ScaleExpr< double, Vector > >
```

The use of templates for encoding such expressions is called the expression template paradigm.
One also refers to *lazy evaluation*, since the evaluation happens later, just when the result is written to the destination vector.



The ASC-bla library implements such expression templates for vectors, in the
[branch *expr*](https://github.com/TUWien-ASC/ASC-bla/tree/expr).

Here we have the base class template `VecExpr` for all vector expressions,
the `SumVecExpr`, and the `operator+` for creating such a SumVecExpr.

```cpp
template <typename T>
  class VecExpr
  {
  public:
    auto derived() const { return static_cast<const T&> (*this); }
    size_t size() const { return derived().size(); }
    auto operator() (size_t i) const { return derived()(i); }
  };

template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SumVecExpr (TA _a, TB _b) : a(_a), b(_b) { }

    auto operator() (size_t i) const { return a(i)+b(i); }
    size_t size() const { return a.size(); }      
  };

template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
{
  return SumVecExpr(a.derived(), b.derived());
}
```

The fancy trick is that `SumVecExpr` derives from the base class `VecExpr`, and gives itself as the template argument to the base. So, we can statically up-cast the base to the derived class.
This idiom is known as
[curiously recurring template pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) or 
[Barton–Nackman trick](https://en.wikipedia.org/wiki/Barton–Nackman_trick#:~:text=The%20idiom%20is%20characterized%20by,recurring%20template%20pattern%20(CRTP).&text=The%20Barton–Nackman%20trick%2C%20then,to%20deal%20with%20such%20ambiguities).
In the breaking work by Todd Veldhuizen [Expression Templates](./Veldhuizen.pdf) the expression templates paradigm for vector operations was introduced. However, back in 1995, it was too much for compiler technology.


If we call the call operator `operator()(size_t)` of an `VecExpr<T>` object, it downcasts to `T`, and calls the call operator there. In this example the `operator()` of a `SumVecExpr` calls the `operator()` of both of its members.

How can this `operator+` be applied to `Vector`s ? Do we have to define all combinations of `operator+([Vector|VecExpr], [Vector|VecExpr])` ? We can avoid it by letting `Vector` derive from `VecExpr`. However, we don't want to copy the vector into the `SumVecExpr`. We could do it by using references - or, alternatively, we introduce a **view** of a vector, a `VectorView`.

## VectorView

Consider the following class hierarchy:

```cpp
template <typename T>
class VecExpr;

template <typename T>
class VectorView : public VecExpr<VectorView<T>>
{
protected:
  size_t m_size;
  T * m_data;
public:
  VectorView (size_t size, T * data)
    : m_size(size), m_data(data) { }

T & operator()(size_t i) { return m_data[i]; }
};

template <typename T>
class Vector : public VectorView<T>
{
 public:
   Vector (size_t size)
     : VectorView (size, new T[size]) { }
   ~Vector() { delete [] m_data; }
};

```

The `VectorView` is a slim class, which can be easily copied using the default copy constructor just copying size and the data pointer. Thus, a `VectorView` can be used as a call-by-value function argument (rather than a reference). All memory allocation/freeing happens in the derived class `Vector`.


A `VectorView` allows also to access a range of a vector:
```cpp
class VectorView {
  ...
  VectorView range(size_t first, size_t next)
    { return VectorView(next-first, m_data+first); }
}
```
With this we can, for example, zero elements with indices in semi-open interval $[10,15)$ via `vec.Range(10,15) = 0`.

A generalization of a `VectorView` allows vectors whose value don't have to lie consecutively in memory. For that we provide the `dist` member variable:
```cpp
class VectorView {
  size_t m_size;
  size_t m_dist;
  T * m_data;
  ... 
  T & operator()(size_t i) { return m_data[i*m_dist]; }  
}
```

Ok, this is more general - but the index calculation comes with some price, which we do not want to pay if we do not need it. As a solution we define the `dist_` variable of a template type, which is set to `std::integral_constant<1>` as a default. The compiler can easily optimize out the multiplication with a constant 1:

```cpp
template <typename T, typename TDIST=std::integral_constant<1>>
class VectorView {
  size_t m_size;
  TDIST m_dist;
  T * m_data;
  ... 
  T & operator()(size_t i) { return m_data[i*m_dist]; }  
}
```

## Exercise

  * Merge the expr - branch from TUWien-ASC/ASC-bla into your main branch [instructions](inst_merge.md)

  * Implement a `MatrixView` class

    ```cpp
    template <typename T, template ORDERING>
    class MatrixView {
      size_t m_rows, m_cols, m_dist;
      T * m_data;
    }
    ```
    Index calculation is `i*m_dist+j` in the row-major case, and `i+j*m_dist` in the col-major case.

    Let your `Matrix` class derive from `MatrixView`. Initialize the `m_dist` variable either with the width, or the height of the matrix, depending on row or column-major storage.

  * Introduce expression templates for matrices, including
    - MatExpr + MatExpr -> MatExpr
    - MatExpr * MatExpr -> MatExpr
    - MatExpr * VecExpr -> VecExpr

    To get the type of the summation variable in the mat-mat or mat-vec products you may use

        using elemtypeA = std::invoke_result<TA,size_t,size_t>::type;
        using elemtypeB = std::invoke_result<TB,size_t>::type;
        using TSUM = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());
        TSUM sum = 0;

  * Implement `matrix.row(i)` and `matrix.col(j)` methods returning a `VectorView` of individual rows and columns.

  * Implement `matrix.rows(first,next)` and `matrix.cols(first,next)` methods returning a `MatrixView` of a range of rows/cols.

  * Implement a `transpose(matrix)` function resulting in a `MatrixView` of opposite ordering
  
  * Simplify the `inverse` function using these new features
  
## Pitfalls


  * What happens with the statement `A = A*B`, where A and B are two matrices ?

  * What is the computational complexity for `A = B*C*D`, where `A`, `B`, `C`, `D` are $n \times n$ matrices ? 

  * what happens here ? 
  ```cpp
  auto func() {
    Matrix a(10,10);
    return a.row(3); 
  }
  ```

    
## How good is it ? 

But can the compiler really generate good code from all of these nested functions and expression objects ? Yes ! It is important that the compiler can inline all the functions, sees the whole flow of data, and optimizes everything as a single function.

To verify what the compiler generates, we can have a look into the generated assembly code. There is an online tool [Compiler Explorer](https://godbolt.org/z/qePEhvaov). You copy in the source code, and it immediately displays the generated assembly code. It allows choosing between a lot of compilers, versions and provided flags.

If you scroll down within the left window you find two functions `MyFunc`, and `MyFunc2`. One uses expression templates, the other one handwritten C-code. In the right window you see the generated assembly code. You can identify a loop (with compare `cmp` and not-equal branching `jne`). You find one addition `addsd` and one multiplication `mulsd` within the loop. You see that the two generated codes are identical, there is no overhead coming from the expression templates.



