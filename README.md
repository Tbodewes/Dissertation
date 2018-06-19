# Dissertation
R code for dissertation on Bayesian network structure learning under missing data using node-average likelihood

Key files:
- Scoring.R: implements node-average likelihood (see Balov, 2013) for conditional Gaussian networks (CGN) and related score functions
- Fit DAG.R: implements local search algorithms for determining a locally optimal CGN for given data and an exact algorithm for when the node order is known and the maximum number of parents is restricted (unsuited for large networks)
- Tests.R: implements tests for scoring and local search

Other files (not up to date with changes in other code):
- Balov alarm experiment.R: implements the computational experiment for the alarm network from Balov (2013)
- 2-node experiment.R: implements two node experiment from Balov (2013)
