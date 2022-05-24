# Implementation for the Stochastic Actor-oriented Models

The Stochastic Actor-Oriented Model (SAOM) for network change (Snijders, 1996) is a model for the dynamical change in networks over time.

We fit three main models to our two networks (arms trade and conventional trade) here:

1. A SAOM of the two networks without interactions between the networks (simple)
2. A model like 1. but we additionally allow for between network interactions (multilevel)
3. A variation of 2. that apply the same specification to sliding windows of 4 years over the whole observation period (sw)

The fitting procedures are in the scripts `fit_saom_*.R`

Results and goodness-of-fits are computed in `results_saom_*.R`
