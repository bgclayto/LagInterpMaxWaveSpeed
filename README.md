# Max Wave Speed Speed Estimator for an Arbitrary EOS in the Lagrangian Frame

## Summary

This program is used to compute or estimate a maximum wave speed for
the Riemann problem of the compressible Euler equations with an
arbitrary equation of state (EOS). This estimation is a slightly more generalized version of the one presented in:

[B. Clayton, J.-L. Guermond, B.
Popov, 2022](https://epubs.siam.org/doi/abs/10.1137/21M1414097). *Invariant
Domain-Preserving Approximations for the Euler Equations with
Tabulated Equation of State*. SIAM Journal on Scientific Computing. 44(1), A444-A470.

and is concerned with Lagrangian frame rather than the Eulerian frame. In particular, the program computes
   
$$
\widehat{\lambda}_{\max} := \max_{Z \in \{L, R\}} \{ \frac{a_Z}{\tau_Z} \sqrt{1 + \frac{\gamma_Z + 1}{2\gamma_Z} \Big( \frac{\widehat{p}^{\*} - p_Z}{p_Z + p_\infty}, 0 \Big)} \}
$$

where $\widehat{p}^*$ is either an upper bound on $p^*$ or is identically $p*$ up to some tolerance ($10^{-9}$ for example).  

## Compile 

In order to compile the program, you need to have `cmake` (at least version
3.5) and `make` installed. 

In `LagInterpMaxWaveSpeed` run the following commands:  
```
mkdir build
cd  build
make release
```
Alternatively, you can run `make debug` for more error checking flags. 

## Execute

Afterwards, To execute the program, run the following commands
```
cd ..
./lambda [argument]
```
where `[argument]` can be replaced by `true` which uses the iterative
method to compute $\lambda_{\text{max}}$ or `false` to use the upper
estimate. If no argument is provided, then the devault value is set to
`true`, i.e. the iterative method is used.

## Test Problems

When `lambda` is run, the max wave speed computation is executed for a
number of test problems given in `data`. These can then be compared
with the true values stored in `output_ref`. Note that the cases 1 -
14 use the ideal gas law while cases 15 - 18 use the van der Waals
EOS.

## Modifications and Details

Note that interpolation parameters `b_covolume`, `q`, and `p_infty`
are assigned to be `0.d0` in `arbitrary_eos.f90`. So, if we would like
to interpolate when the oracle is given by the van der Waals EOS, then
we should assign `b_covolume` to be the same value provided by van der
Waals EOS, `b_vdw`. However, for simplicity, all computations were run
with all interpolation parameters set to zero.
