# Hidden Markov Models and Their Application in Modeling Rainfall Occurrence

MSc thesis in Statistics, Athens University of Economics and Business (AUEB), March 2025.
**Author:** Konstantinos Grammenos · **Supervisor:** Prof. Panagiotis Besbeas

📄 [Full thesis (PDF)](./Grammenos_2025.pdf)

## Overview

This thesis uses Hidden Markov Models (HMMs) to model daily rainfall occurrence over a
90-day wet season (Feb–Apr), across 24 years (1975–2002), at 10 weather stations in the
state of Ceará, Northeast Brazil.

The analysis progresses through several modeling stages:

1. **Introductory binary HMM** on the Old Faithful Geyser eruption data, as a worked
   example of state-dependent binary time series.
2. **Homogeneous HMMs** (2–5 states) fitted to the multivariate rainfall data, identifying
   4 hidden states: a wet state, a dry state, and two "transitional" states with contrasting
   north/south spatial rainfall patterns.
3. **Non-Homogeneous HMM (NHMM)**, incorporating a GCM-simulated seasonal rainfall
   anomaly covariate into the transition probability matrix.
4. **Covariate-in-observation models**, where the same climate covariate is added to the
   emission (Bernoulli) probabilities via logistic regression, with both shared and
   station-specific slope specifications.
5. **Independent seasonal sequences**, modeling each of the 24 years as its own sequence
   rather than one continuous 2160-day chain.

Model selection throughout is based on AIC/BIC and (for the baseline HMM) cross-validated
log-likelihood.

## Repository contents

| File | Description |
|---|---|
| `Grammenos_2025.pdf` | Full thesis text |
| `precip_forecast_code.R` | R implementation of a Non-Homogeneous Poisson HMM (EM algorithm with conjugate-gradient M-step) used for the precipitation-anomaly-driven modeling and forecasting |

## Data

Rainfall occurrence data for the 10 Ceará stations (1975–2002) were originally compiled by
Robertson, Kirshner & Smyth (2003), sourced from FUNCEME (Fundação Cearense de
Meteorologia e Recursos Hídricos). The GCM covariate is the ECHAM 4.5 simulated seasonal
mean rainfall anomaly (Roeckner et al., 1996). Data are not included in this repository —
contact the original sources for access.

## Key findings

- 4 hidden states best describe rainfall dynamics: wet, dry, and two north/south contrasting states.
- The GCM seasonal rainfall anomaly significantly affects transition dynamics between states.
- Incorporating the covariate into station-level rainfall probabilities (rather than just
  transitions) shows a consistent, near-uniform effect across stations.
- Modeling each season as an independent sequence (rather than one continuous chain)
  gives a marginally better fit, consistent with year-to-year shifts in ITCZ/SAMS behavior.

## References

See the bibliography in the thesis PDF, in particular Zucchini, MacDonald & Langrock (2017),
*Hidden Markov Models for Time Series*, and Hughes & Guttorp (1994) on non-homogeneous HMMs for precipitation.
