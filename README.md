# Dissertation
Bayesian network structure learning from incomplete data using Node-Average Likelihood
R code and experimental results

R code:
- Scoring.R: Implements NAL scoring.
- Fit DAG.R: Implements algorithms for fitting with known and unknown ordering. Also implements parametric and structural EM and a function for bootstrapping DAGs.
- Experiments.R: Implements functions necessary to run computational experiments.
- Analysis.R: Code to configures and run experiments and process results.
- MEHRA functions.R: Implements functions necessary for MEHRA analysis.
- MEHRA analysis.R: Script for MEHRA analysis.
- Tests.R: Integration tests for key functions.

Experimental results:
- Result_ordered.RDS: results for experiment with known node ordering
- Result_unordered.RDS: results for experiment with unknown node ordering
- Result_em.RDS: results for comparison of NAL and Structural EM
