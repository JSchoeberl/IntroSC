## Refresh our C++ knowledge

We assume you have a basic knowledge of C++, and repeat only shortly some features of the C++ programming language. Here is a list of [C++ books](https://stackoverflow.com/questions/388242/the-definitive-c-book-guide-and-list).




```cpp
template <typename T>
class Vector
{
  size_t size_;
  T * data_;

public:
  Vector (size_t size) 
    : size_(size), data_(new T[size]) { ; }

  Vector (const Vector & v)
    : Vector(v.Size())
  {
    *this = v;
  }

  Vector (Vector && v)
    : size_{0}, data_(nullptr)
  {
    std::swap(size_, v.size_);
    std::swap(data_, v.data_);
  }

  ~Vector () { delete [] data_; }

  Vector & operator= (const Vector & v2)
  {
    for (size_t i = 0; i < size_; i++)
      data_[i] = v2(i);
    return *this;
  }

  Vector & operator= (Vector && v2)
  {
    for (size_t i = 0; i < size_; i++)
      data_[i] = v2(i);
    return *this;
  }

  size_t Size() const { return size_; }
  T & operator() (size_t i) { return data_[i]; }
  const T & operator() (size_t i) const { return data_[i]; }
};
```

* We define a class template `Vector`. The type of the vector elements is given by template argument `T`. For example, one can define a vector of double precision entries as `Vector<double>`.

* A Vector object can be created by one of its constructors, and is destroyed by the destructor.

* We use the *move-semantics* introduced with C++11. If a Vector is initialized from a return-value Vector (which will be destroyed anyway), the move-constructor can steal resources from the temporary vector.

* All data members are private, what means invisible form outside. The interaction with an Vector object happens only via its methods (`Size`, assignment `operator=` and call `operator()`). We have chosen () instead of [] for consistency with the `Matrix` class.



The function `operator+` for adding two vectors of type `T` is as follows:

```cpp
template <typename T>
Vector<T> operator+ (const Vector<T> & a, const Vector<T> & b)
{
  Vector<T> sum(a.Size());
  for (size_t i = 0; i < a.Size(); i++)
    sum(i) = a(i)+b(i);
  return sum;
}  
```

