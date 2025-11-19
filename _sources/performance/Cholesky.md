# Efficient Cholesky factorization

The Cholesky factorization of a symmetric, positive definite matrix $A$ is

$$
A = L L^T,
$$

where $L$ is a lower triangular matrix. After computing the factorization,
the solution of a linear system $A x = b$ can be found by forward/backward substitution.

Very similar is the `LDL` factorization

$$
A = L D L^T
$$

where $D$ is a diagonal matrix, and $L$ is a lower, normalized (i.e. $L_{ii}=1$) matrix.
The advantage is to avoid computing square roots, and it (may) also work for non positive definite matrices.

A simple algorithm for computing the LDL-factorization is

```cpp
void CalcLDL (MatrixView<T,ORD> A)
{
  size_t n = A.rows();
  for (size_t i = 0; i < n-1; i++) {
     T inv_dii = 1.0 / A(i,i);
     for (size_t j = i+1; j < n; j++) {
       T hji = A(j,i)
       A(j,i) = hji*inv_dii:
       for (size_t k = i+1; k <= j; k++)
         A(j,k) -= hji * A(k,i)
  }
}
```

If we want to optimize that algorithm, we have to redo the optimization of micro-kernels, and the cache optimization.


## A recursive computation of the Cholesky factorization

The Cholesky-factorization of a $1 \times 1$ matrix is just taking the square-root of $a_{11}$.
For matrix sizes $n \geq 2$ we split the matrix into $2 \times 2$ block matrices:


\begin{eqnarray*}
\left( \begin{array}{cc}
   A_{11} & * \\
   A_{21} & A_{22}
      \end{array} \right) & = &

\left( \begin{array}{cc}
  L_{11} & 0 \\
  L_{21} & L_{22}
    \end{array} \right)
\left( \begin{array}{cc}
  L_{11}^T & L_{21}^T \\
   0  & L_{22}^T
    \end{array} \right) \\
& = & 
\left( \begin{array}{cc}
  L_{11} L_{11}^T & * \\
  L_{21} L_{11}^T &  L_{21} L_{21}^T + S
    \end{array} \right)
\end{eqnarray*}
with $S = L_{22} L_{22}^T$. Here, $*$ means defined by symmetry.


The block-Cholesky factorization can be computed by 4 steps:

1. Compute Cholesky factoriaztion for $A_{11} = L_{11} L_{11}^T$
2. Compute $L_{21}$ by solving a triangular system $L_{21} L_{11}^T = A_{21}$
3. Compute $S = A_{22} - L_{21} L_{21}^T$
4. Compute Cholesky factorization of $S = L_{22} L_{22}^T$.

Step 1 and step 4 are Cholesky factorizations of matrices of size $n/2$.
Step 3 is a matrix-matrix product, which we can compute with optimal performance.

Step 2 is solving a triangular system for many vectors $(x_1, \ldots x_n)$ stored in a matrix $X$:

$$
L X = Y
$$

Again, we write $L$ as a $2 \times 2$ block matrix, and split $X$ and $Y$ vertically:

$$
\left( \begin{array}{cc}
  L_{11} & \\
  L_{21} & L_{22} 
  \end{array} \right)
\left( \begin{array}{cc}
  X_1 \\ X_2
  \end{array} \right)
= 
\left( \begin{array}{cc}
  Y_1 \\ Y_2
  \end{array} \right)
$$

This is done in 3 steps:
1. Compute $X_1$ by solving $L_{11} X_1 = Y_1$
2. Compute $\widetilde Y_2 = Y_2 - L_{21} X_1$
3. Compute $X_2$ by solving $L_{22} X_2 = \widetilde Y_2$.

Steps 1 and 3 are problems of the same type, but for half of the problem size.
Note that the right hand side matrix $Y$ can be overwritten by the solution $X$, no additional memory needs to be allocated.


```cpp
// compute X <-- Inv(L) * X
// L ... n \times n
// X ... n \times m
template <typename TL, Ordering OL, typename TX, Ordering OX>
void SolveTriangular (MatrixView<TL,OL> L, MatrixView<TX,OX> X) {
  size_t n = L.rows();
  if (n == 0) return;
  if (n == 1) {
    X *= 1/L(0,0);
    return;
  }

  size_t n2 = n/2;
  auto L11 = L.rows(0,n2).cols(0,n2);
  auto L21 = L.rows(n2,n).cols(0,n2);
  auto L22 = L.rows(n2,n).cols(n,n);

  auto X1 = X.rows(0,n2);
  auto X2 = X.rows(n2,n);

  SolveTriangular (L11, X1);
  X2 -= L21 * X1;
  SolveTriangular (L22, X2);
}
```






The `LDL` is a little bit more tricky, since 

\begin{eqnarray*}
\left( \begin{array}{cc}
   A_{11} & * \\
   A_{21} & A_{22}
      \end{array} \right) & = &
\left( \begin{array}{cc}
  L_{11} & 0 \\
  L_{21} & I
    \end{array} \right)
\left( \begin{array}{cc}
  D_{11} & 0 \\
  0 & S 
    \end{array} \right)
\left( \begin{array}{cc}
  L_{11}^T & L_{21}^T \\
  0 & I
    \end{array} \right) \\
& = & 
\left( \begin{array}{cc}
  L_{11} D_{11} L_{11}^T & * \\
  L_{21} D_{11} L_{11}^T &  L_{21} D_{11} L_{21}^T + S
    \end{array} \right)
\end{eqnarray*}

