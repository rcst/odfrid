# odfrid - *o*rigin-*d*estination estimation *f*rom *rid*ership data

## Notes

## TODO

- [ ] Conditional Sampling of OD vector: Metropolis-Hasting
- [x] implement generation of routing matrix
- [ ] implement solution-check function `H` (= inv(A) * x)
- [ ] implement softmax transformation for alighting probabilities vector (function of \rho and G_i)
- [ ] implement rank-check function for matrix `G` (assumption is that `G` has a low-rank structure)
- [ ] fake-data generation for testing
- [ ] import sample data or implement
- [ ] write function to plot bus passenger load plots from OD vector

### Conditional Sampling of OD vector: Metropolis-Hasting

NOTE: 
Sampling from current parameter settings $\Theta$ is achieved by means of the 
acceptance probabilities dependence on the model parameters

[x] (cpp) calculate vector $w$ from stacked vectors $u$ and $v$ iteratively 
$w_i = w_{i-1} + u_i - v_i$, with the initial condition $w_0 = 0$
[ ] (cpp) implement function that calculates a valid but random OD-vector $y$
    - package the output with vector $z$ that is calculated simulatenously to avoid
    recalculation of
    - package the output of the Markov chain transition probabilities $\pi$ with
      it
