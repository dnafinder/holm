# holm
Holm-Sidak t-test: a routine for multiple t-test comparisons </br>
This file is applicable for equal or unequal sample sizes

Syntax: 	HOLM(X,GROUP,CTRL,ALPHA,TAIL)
     
    Inputs:
          X: data vector
          GROUP - specifies grouping variables G. Grouping variables must
          have one column per element of X.
          CTRL: The first sample is a control group (1); there is not a
          control group (0). (default=0).
          ALPHA: significance level (default = 0.05).
          TAIL: 2 (2-tailed, default);
               -1 (left tailed)
                1 (right tailed)
    Outputs:
          - Mean and Standard Deviation vectors.
          - degrees of freedom and combined variance.
          - p-value for each comparison.
          - alpha value corrected by Sidak procedure.
          - whether or not Ho is rejected.

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2006). Holm-Sidak t-test: a routine for multiple t-test comparisons.
http://www.mathworks.com/matlabcentral/fileexchange/12786
