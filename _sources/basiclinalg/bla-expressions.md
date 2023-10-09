## Expression templates

Think of a code like



```cpp
z = x + 3*y;
```

where `x,y,z` are vectors. What happens ? First, the `operator* (double, Vector)` is called, it generates a temporary vector for 3*y, and then the `operator+ (Vector, Vector)` is called to compute a new vector for the sum. Finally, the later vector is copied into the memory of `z`.

Writing an old-style C function 
```cpp
for (size_t i = 0; i < z.Size(); i++)
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



The ASC-bla library implements such expression templates for vectors. Git-clone from
[branch *expr*](https://github.com/TUWien-ASC/ASC-bla/tree/expr).

Here we have the base class template `VecExpr` for all vector expressions,
the `SumVecExpr`, and the `operator+` for creating such a SumVecExpr.

```cpp
template <typename T>
  class VecExpr
  {
  public:
    auto View() const { return static_cast<const T&> (*this); }
    size_t Size() const { return View().Size(); }
    auto operator() (size_t i) const { return View()(i); }
  };

template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a_;
    TB b_;
  public:
    SumVecExpr (TA a, TB b) : a_(a), b_(b) { }
    auto View () { return SumVecExpr(a_, b_); }

    auto operator() (size_t i) const { return a_(i)+b_(i); }
    size_t Size() const { return a_.Size(); }      
  };

template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
{
  return SumVecExpr(a.View(), b.View());
}
```

The fancy trick is that `SumVecExpr` derives from the base class `VecExpr`, and gives itself as the template argument to the base. So, we can statically up-cast the base to the derived class.
This idiom is known as
[curiously recurring template pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) or 
[Barton–Nackman trick](https://en.wikipedia.org/wiki/Barton–Nackman_trick#:~:text=The%20idiom%20is%20characterized%20by,recurring%20template%20pattern%20(CRTP).&text=The%20Barton–Nackman%20trick%2C%20then,to%20deal%20with%20such%20ambiguities).
In the breaking work by Todd Veldhuizen [Expression Templates](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ca2f8a9b7407de039957a358f995265ec8b769a9) the expression templates paradigm for vector operations was introduced. However, back in 1995, it was too much for compiler technology.


If we call the call operator `operator()(size_t)` of an `VecExpr<T>` object, it upcasts to `T`, and calls the call operator there. In this example the `operator()` of a `SumVecExpr` calls the `operator()` of both of its members.

How can this `operator+` be applied to a `Vector` ? Do we have to define call combinations of `operator+([Vector|VecExpr], [Vector|VecExpr])` ? We can avoid it by letting `Vector` derive from `VecExpr`. However, we don't want to copy the vector into the `SumVecExpr`. We could do it by using references - or, alternatively, we introduce a **view** of a vector, a `VectorView`.

### VectorView

Consider the following class hierarchy:

```cpp
class VecExpr;

template <typename T>
class VectorView : public VecExpr<VectorView<T>>
{
protected:
  size_t size_;
  T * data_;
public:
  VectorView (size_t size, T * data)
    : size_(size), data_(data) { }
  VectorView (const Vec
  ~VectorView() { }
  T & operator()(size_t i) { return data[i]; }
};

template <typename T>
class Vector : public VectorView<T>
{
 public:
   Vector (size_t size)
     : VectorView (size, new T[size]) { }
   ~Vector() { delete [] data_; }
};

```

The `VectorView` is a slim class, which can be easily copied using the default copy constructor just copying size and the data pointer. Thus, a `VectorView` can be used as a call-by-value function argument (rather than a reference). All memory allocation/freeing happens in the derived class `Vector`.


A `VectorView` allows also to access a range of a vector:
```cpp
class VectorView {
  ...
  VectorView Range(size_t first, size_t next)
    { return VectorView(next-first, data+fisrt); }
}
```
With this we can, for example zero elements  with indices $10 \leq i < 20$ via `vec.Range(10,20) = 0`.

A generalization of a `VectorView` is to allow vectors whose value don't have to lie consecutatively in memory. For that we provide the `dist` member variable:
```cpp
class VectorView {
  size_t size_;
  size_t dist_;
  T * data_;
  ... 
  T & operator()(size_t i) { return data[i*dist_]; }  
}
```

Ok, this is more general - but the index calculation comes with some price, which we do not want to pay if we do not need it. As a solution we define the `dist_` variable from a template type, which is set to `std::integral_constant<1>` as a default. The compiler can easily optimize out the multiplication with a constant 1:

```cpp
template <typename T, typename TDIST=std::integral_constant<1>>
class VectorView {
  size_t size_;
  TDIST dist_;
  T * data_;
  ... 
  T & operator()(size_t i) { return data[i*dist_]; }  
}
```

### Excercise

  * Implement a `MatrixView` class

    ```cpp
    template <typename T, template ORDERING>
    class MatrixView {
      size_t height, width, dist;
      T * data_;
    }
    ```
    Index calculation is `i*dist+j` in the row-major case, and `i+j*dist` in the col-major case.

  * Introduce an expression-template hierarchy for matrices, including
    - Matrix + Matrix -> MatExpr
    - Matrix * Matrix -> MatExpr
    - Matrix * Vector -> VecExpr


  * Implement `matrix.Row(i)` and `matrix.Col(j)` methods returning a `VectorView` of individual rows and columns.

  * Implement a `Transpose(matrix)` function resulting in a `MatrixView` of opposite ordering
  
   
### How good is it ? 

But can the compiler really generate good code from all of this nested functions and expression objects ? Yes ! It is important that the compiler can inline all the functions, sees the whole flow of data, and optimizes everything as a single function.

To verify what the compiler generates, we can have a look into the generated assembly code. There is an online tool [Compiler Explorer](https://godbolt.org/z/qePEhvaov). You copy in the source code, and it immediately displays the generated assembly code. It allows to choose between a lot of compilers, versions and provided flags.

If you scroll down within the left window you find two functions `MyFunc`, and `MyFunc2`. One uses expresion templates, the other one hand-written C-code. In the right window you see the generated assemply code. You can identify a loop (with compare `cmp` and not-equal branching `jne`. You find one addition `addsd` and one multiplication `mulsd` within the loop. You see that the two generated codes are identic, there is no overhead coming from the expression templates.



