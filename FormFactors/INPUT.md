# WIMpy_NREFT


## NREFT form factors

The nuclear form factors are written in the form:
$$
F_X^{(N, N')} = \mathrm{e}^{-2y}\sum_{k} c_k y^k\,.
$$

The definition of the dimensionless momentum variable y can be found in e.g. Appendix A.3 of [arXiv:1203.3542](https://arxiv.org/abs/1203.3542).

The form factors tables give the coefficients $c_k$ from $k = 0$ to $k = 7$ on each row. For each form factor $F_X$ there are 4 rows, corresponding to:

$F_X^{(p,p)}\,,$

$F_X^{(p,n)}\,,$

$F_X^{(n,p)}\,,$

$F_X^{(n,n)}\,.$

The different form factors are listed in the following order (with 4 rows for each, as described above):
$$
F_M, \,\,F_{\Sigma'}, \,\, F_{\Sigma''}^{(p,p)}, \,\, F_{\Delta}, \,\, F_{\Phi''}, \,\, F_{M,\Phi''}, \,\, F_{\Sigma',\Delta}\,.
$$

Note that the form factor $F_{\tilde{\Phi}'}$ is not (yet) included, and therefore the results are only valid for operators up to $\mathcal{O}_{11}$.

### References

NREFT form factors for Xe, F and I are taken from [arXiv:1203.3542](https://arxiv.org/abs/1203.3542). The form factor for C is taken from [arXiv:1501.03729](https://arxiv.org/abs/1501.03729).

Note that the $W$ functions of arXiv:1501.03729 (and others) are related to the form factors $F$ by:

$$
F_X^{(N,N')} = \frac{4\pi}{2 J + 1} W_X^{(N, N')}\,,
$$
where J is the nuclear spin. Check out Eq. 76 of arXiv:1203.3542 for how to convert between the nucleon $(N, N')$ and isospin $(\tau, \tau')$ bases.