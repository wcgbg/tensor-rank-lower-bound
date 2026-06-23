# References used for `prev_results.csv`

## External bibliographic and database sources

1. C. M. Fiduccia and Y. Zalcstein, “Algebras Having Linear Multiplicative Complexities,” Journal of the ACM, 1977. DOI: 10.1145/322003.322014.
2. N. H. Bshouty and M. Kaminski, “Multiplication of Polynomials over Finite Fields,” SIAM Journal on Computing 19(3), 1990. DOI: 10.1137/0219029.
3. S. Winograd, “Some Bilinear Forms Whose Multiplicative Complexity Depends on the Field of Constants,” Mathematical Systems Theory 10, 1976/77. DOI: 10.1007/BF01683270.
4. S. Ballet, J. Chaumine, J. Pieltant, M. Rambaud, H. Randriambololona, and R. Rolland, “On the Tensor Rank of Multiplication in Finite Extensions of Finite Fields and Related Issues in Algebraic Geometry,” Russian Mathematical Surveys 76(1), 2021. DOI: 10.1070/RM9928. Preprint: arXiv:1906.07456.
5. S. Ballet, A. Bonnecaze, and B. Pacifico, “Multiplication in Finite Fields with Chudnovsky-Type Algorithms on the Projective Line,” arXiv:2007.16082v2, 2021.
6. M. Grassl, “Bounds on the Minimum Distance of Linear Codes and Quantum Codes,” CodeTables.de, online database. https://www.codetables.de/.
7. A. Alder and V. Strassen, “On the Algorithmic Complexity of Associative Algebras,” Theoretical Computer Science 15, 201–211, 1981.
8. M. Bläser, “A Complete Characterization of the Algebras of Minimal Bilinear Complexity,” SIAM Journal on Computing 34(2), 277–298, 2004. DOI: 10.1137/S0097539703438277.
9. A. Averbuch, Z. Galil, and S. Winograd, “Classification of All the Minimal Bilinear Algorithms for Computing the Coefficients of the Product of Two Polynomials Modulo a Polynomial. Part II: The Algebra G[u]/<u^n>,” Theoretical Computer Science 86(2), 143–203, 1991. DOI: 10.1016/0304-3975(91)90017-V.
10. L. Wang and Q. Wang, “On Explicit Factors of Cyclotomic Polynomials over Finite Fields,” Designs, Codes and Cryptography 63, 87–104, 2012; preprint arXiv:1011.4857.
11. A. Boripan and S. Jitman, “Revisiting the Factorization of x^n+1 over Finite Fields,” arXiv:2009.09601, 2020.

## Derived ingredients recorded in the CSV

- Griesmer bound: used as a derived coding-theoretic lower-bound calculation; the CSV does not encode a separate bibliographic citation for it.
- Factor-field quotient monotonicity: used as a derived cyclic-convolution lower-bound transfer.
- CRT decomposition into extension-field components: used as a derived cyclic-convolution upper-bound construction.
- Linear folding modulo x^n-1 of a full polynomial-product algorithm: used as a derived cyclic-convolution upper bound.
- Karatsuba recursion after zero-padding: used as a derived full-product upper-bound construction.
- Trivial scalar multiplication: the n=1 base case.
- Exact truncated-product rank: F_q[x]/(x^N) is a one-factor local quotient; the table records R=2N-1 for all listed q,N.
- Cyclotomic-order factor degrees: for squarefree parts, Φ_r over F_q splits into φ(r)/ord_r(q) irreducibles of degree ord_r(q); repeated p-powers are tracked as multiplicities.
- Negacyclic characteristic-2 identity: x^N+1=x^N-1, so the negacyclic row is copied from the improved cyclic row.
- Local-extension component upper bound: a factor f(x)^e with deg(f)=d is bounded by multiplying in F_{q^d}[t]/(t^e), costing (2e-1) times the chosen upper bound for multiplication in F_{q^d}; component full-product folding is also compared.
- Cyclic repeated-factor improvement: after adding truncated-product exact ranks, repeated-root cyclic cases are updated by CRT/local-component upper bounds when they improve the earlier folding upper bound.

## Reference-cell registry

The generator script contains the exact `lb_ref` and `ub_ref` text emitted for each row.
