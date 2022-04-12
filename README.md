# Bayesian-sparse-PCE

This repository contains supporting material for the paper entitled *The sparse
Polynomial Chaos expansion: a fully Bayesian approach with joint priors on the
coefficients and global selection of terms* by Paul-Christian Bürkner, Ilja
Kröker, Sergey Oladyshkin, and Wolfgang Nowak. The CO2 data set, which is used
in one of our case studies, has been already published elsewhere (Köppel et al.,
2017) and should be cited as such in case of re-use in accordance with its
license.

In the code, we are using partially different notation than in the paper,
as a result of how this project has evolved over time. Below, we list the
main concepts and their corresponding notation in paper and code:

- Number of input variables: *N* in the paper, `M` in the code
- Number of training points: *T* in the paper, `N` or `mc` in the code
- Maximal polynomial degree: *d* in the paper, `p` in the code
- Total number of polynomials: *M* in the paper, `P` in the code
- Input variables: <img src="https://render.githubusercontent.com/render/math?math=\omega"> in the paper, `x` in the code
- Evaluated polynomial matrix: <img src="https://render.githubusercontent.com/render/math?math=\Psi"> in the paper, `X` in the code
- Polynomials of the PCE: *c* in the paper, `b` in the code


## References

Köppel, M., Franzelin, F., Kröker, I., Oladyshkin, S., Wittwar, D., Santin, G.,
Barth, A., Haasdonk, B., Nowak, W., Pflüger, D., and Rohde, C. Datasets and
executables of data-driven uncertainty quantification benchmark in carbon dioxide storage. (2017) https://doi.org/10.5281/zenodo.933827
