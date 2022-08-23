<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Make
grid](https://github.com/ecosystem-state/rreas-auto/workflows/R-make-grid/badge.svg)](https://github.com/ecosystem-state/rreas-auto/actions)

[![Estimate
densities](https://github.com/ecosystem-state/rreas-auto/workflows/R-make-indices/badge.svg)](https://github.com/ecosystem-state/rreas-auto/actions)  
[![Run DFA
models](https://github.com/ecosystem-state/rreas-auto/workflows/R-run-dfa/badge.svg)](https://github.com/ecosystem-state/rreas-auto/actions)
<!-- badges: end -->

## Overview

This repository is a demonstration of automatic index generation using
data from RREAS and ERDDAP. The index is generated by applying Dynamic
Factor Analysis (DFA) to the top \~ 20+ species, using samples collected
1990 - present.

## Results

For illustrative purposes we’re only showing the model with 1-trend.

<div class="figure">

<img src="figures/trends.jpeg" alt="Estimated trends for the RREAS community" width="480" />
<p class="caption">
Estimated trends for the RREAS community
</p>

</div>

<div class="figure">

<img src="figures/loadings.jpeg" alt="Estimated loadings for the RREAS community" width="528" />
<p class="caption">
Estimated loadings for the RREAS community
</p>

</div>

<div class="figure">

<img src="figures/fitted.jpeg" alt="Predicted and observed fits to the RREAS data" width="720" />
<p class="caption">
Predicted and observed fits to the RREAS data
</p>

</div>

## bayesdfa

For more on the approach used, check out the [bayesdfa R
package](https://fate-ewi.github.io/bayesdfa/)
