# Expression templates

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
One also referes to *lazy evaluation*, since the evaluation happens later, just when the result is written to the destination vector.



The ASC-bla library implements such expression templates for vectos. Git-clone from
[branch *expr*](https://github.com/JSchoeberl/ASC-bla/tree/expr).

Here we have the base class template `Expr` for all vector and matrix expressions,
the `SumExpr`, and the `operator+` for creating such a SumExpr.

```cpp
template <typename T>
  class Expr
  {
  public:
    auto View() const { return static_cast<const T&> (*this); }
    size_t Size() const { return View().Size(); }
    auto operator() (size_t i) const { return View()(i); }
  };

template <typename TA, typename TB>
  class SumExpr : public Expr<SumExpr<TA,TB>>
  {
    TA a_;
    TB b_;
  public:
    SumExpr (TA a, TB b) : a_(a), b_(b) { }
    auto View () { return SumExpr(a_, b_); }

    auto operator() (size_t i) const { return a_(i)+b_(i); }
    size_t Size() const { return a_.Size(); }      
  };

template <typename TA, typename TB>
  auto operator+ (const Expr<TA> & a, const Expr<TB> & b)
{
  return SumExpr(a.View(), b.View());
}
```

The fancy trick is that `SumExpr` derives from the base class `Expr`, and gives itself as the template argument to the base. So, we can statically up-cast the base to the derived class.
This idiom is known as
[curiously recurring template pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) or 
[Barton–Nackman trick](https://en.wikipedia.org/wiki/Barton–Nackman_trick#:~:text=The%20idiom%20is%20characterized%20by,recurring%20template%20pattern%20(CRTP).&text=The%20Barton–Nackman%20trick%2C%20then,to%20deal%20with%20such%20ambiguities).


If we call the paranthesis operator `operator()(size_t)` of an `Expr<T>` object, if upcasts to `T`, and calls the paranthesis operator there. In this example the `operator()` of a `SumExpr` calls the `operator()` of both of its members.




But can the compiler really generate good code from all of this nested functions and expression objects ? Yes ! It is important that the compiler can inline all the functions, sees the whole flow of data, and optimizes everything as a single function.

To verify what the compiler generates, we can have a look into the generated assembly code. There is a online tool [Compiler Explorer](https://godbolt.org/z/qePEhvaov). You copy in the source code, and it immediately displays the generated assembly code. It allows to choose between a lot of compilers, versions and provided flags.

What we see is ....
