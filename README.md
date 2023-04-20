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
\widehat{\lambda}_{\max} := \max_{Z \in \{L, R\}} \{ \frac{a_Z}{\tau_Z} \sqrt{1 + \frac{\gamma_Z + 1}{2\gamma_Z} \Big( \frac{\widehat{p}^* - p_Z}{p_Z + p_\infty}, 0 \Big)} \}
$$

where $\widehat{p}^{*}$ is either an upper bound on $p^*$ or is identically $p^*$ up to some tolerance ($10^{-9}$ for example).

## Compile 

In order to compile the you need to have `cmake` (at least version 3.5) and `make` installed. 
In `LagInterpMaxWaveSpeed` run the following commands:  

```
mkdir build
cd  build
make release
```

You can also execute `make debug` for more error checking flags. To execute the program, run the following commands:

```
cd run
./a.exe
```

## Test Problems

The simulation is executed for a number of test problems given in `data` compared with the true values stored in `output_ref`.
