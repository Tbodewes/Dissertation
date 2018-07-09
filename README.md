# Dissertation
R code for dissertation on Bayesian network structure learning under missing data using node-average likelihood

Key files:
- Scoring.R: implements node-average likelihood (see Balov, 2013) for conditional Gaussian networks (CGN) and related score functions
- Fit DAG.R: implements local search algorithms for determining a locally optimal CGN for given data and an exact algorithm for when the node order is known and the maximum number of parents is restricted (unsuited for large networks)
- Experiments.R: functions for computational experiments, analyzing SHD from learned graph to true graph and number of queries used in learning, for ordered fitting, general fitting and structural EM. Runs for a set of graphs, a set of sample sizes and a set of missingness probabilities.
- Tests.R: unit tests
