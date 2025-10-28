# Automatic testing

Usually you implement some algorithms, test them - and they work. But often you modify your code later, and it happes that you break code that worked before. It often takes a long time to find these errors.

For that reason it is a good idea to add a test-suite of well chosen tests that verify the individual units of your library.


We are going to use [Catch2](https://github.com/catchorg/Catch2) for setting up such tests.


Now, merge the branch `catch` into your project. It gives you an extra directory `tests` containing the test suite.

```cpp
TEST_CASE( "bla-tests", "[vector size test]" ) {
  Vector<double> x(5);
  Vector<double> y = x;
  REQUIRE(y.size() == 5);
}

TEST_CASE( "bla-tests", "[vector add and scale test]" ) {
  Vector<double> x(5);
  x = 3;
  Vector<double> y = x + x + (-2) * x + x;
  REQUIRE(y(0) == 3); 
}

TEST_CASE( "bla-tests", "[buggy vector add and scale test]" ) {
  Vector<double> x(5);
  x = 3;
  Vector<double> y = x + x + (-2) * x + x;
  REQUIRE(y(0) == 4);  // wrong test !!
}

```

Compiling and running the tests will produce the following output:


    /Users/joachim/texjs/lva/IntroSC/ASC-bla/tests/tests_vector.cpp:27: FAILED:
       REQUIRE( y(0) == 4 )
    with expansion:
      3.0 == 4

    ===============================================================================
    test cases: 3 | 2 passed | 1 failed
    assertions: 3 | 2 passed | 1 failed


Actually, not the vector class, but the test case is buggy.



In vs-code, you can use the testing symbol to run the tests:
```{image} vscode-testing.png
:width: 40%
:align: center
```


One can set up the GitHub server such that with every commit the test-suite is run automatically, see
[Continuous Integraton](https://github.com/resources/articles/devops/ci-cd).
That should help to detect errors at an early stage.



