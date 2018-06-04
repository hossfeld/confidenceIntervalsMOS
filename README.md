# confidenceIntervalsMOS
Confidence Interval Estimators for MOS Values 

The function computes confidence intervals for QoE data on a discrete 5-point rating scale. In particular, the confidence intervals for the Mean Opinion Scores are computed. 
 
Different approaches are used to derive the confidence intervals for MOS values. The first three approaches 1-3 estimate the confidence interval for mean values. The approaches 4-6 consider the fact that user ratings follow a discrete, truncated random variable for a certain TC. Then, the binomial distribution may be used as an upper bound (for the emerging variance due to the binomial sum variance inequality). Accordingly, CI estimators for binomial proportions are then utilized. Then exact CIs of simultaneous confidence intervals for multinomial proportions (7) as well as a bayesian approach (8) are considered.

1. CI for Mean: Bootstrap algorithm
2. CI for Mean: Normal approximation for the MOS
3. CI for Mean: Student-T approximation for the MOS
4. CI for Bino: Wald interval employing normal approximation
5. CI for Bino: Wilson score interval with continuity correction
6. CI for Bino: Clopper-Pearson (using beta distribution)
7. Simultaneous CI: S. Jin, Computing Exact Confidence Coefficients of Simultaneous  Confidence, Intervals for Multinomial Proportions and their Functions. Department of Statistics, Uppsala University, 2013.
8. Jeffrey's interval: H. Jeffreys, The theory of probability. OUP Oxford, 1998.

Input data: 
    alpha  - significance level, default value is 0.05
    y  - user ratings y of size k x n with k test conditions and n subjects

Output data: The return value is a struct containing the following data.
      str - name of the confidence interval estimators
      shortstr - abbreviation (for plots)
      mos - mean opinion score averaged over the n subjects for all k test conditions
      ciUpper - the upper confidence intervall bound
      ciLower - the lower confidence intervall bound
      ciWidth - the size of the confidence interval; note that the CI is not symmetric around the MOS
