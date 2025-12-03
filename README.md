[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/holm&file=holm.m)

# holm

Holm-Sidak procedure for multiple Student's t-tests between independent groups, implemented in MATLAB.

This function performs the Holm-Sidak stepdown procedure using pooled-variance t-tests to compare multiple independent samples, with an optional control group and structured outputs.

## Syntax

The function uses a single intuitive calling style: one data vector per group, followed by optional Name–Value pairs.

Basic forms:

- holm(X1, X2, ..., XK)
- holm(X1, X2, ..., XK, 'Ctrl', CTRL, 'Alpha', ALPHA, 'Tail', TAIL, 'Display', DISPLAY)
- stats = holm(X1, X2, ..., XK, ...)
- [stats, results] = holm(X1, X2, ..., XK, ...)

where:

- X1, ..., XK are numeric vectors (row or column) containing the observations for each of the K groups.
- CTRL is a logical-like flag indicating whether the first group X1 is to be treated as a control group.
- ALPHA is the global significance level (scalar in (0,1), default 0.05).
- TAIL is the tail of the test:
  - 2  → two-sided (default)
  - -1 → left-tailed
  - 1  → right-tailed
- DISPLAY is a logical-like flag controlling command-window output.

Logical-like values accepted for CTRL and DISPLAY include:

- logical: true, false
- numeric: 1, 0
- char/string (case-insensitive): "on", "off", "yes", "no", "true", "false"

## Inputs

### Required inputs

- X1, ..., XK  
  - Numeric vectors (row or column), real, finite, non-NaN, non-empty.  
  - Each Xi contains the data for the i-th group.  
  - Group sizes may be unequal.  
  - At least two groups are required (K ≥ 2).

### Name–Value pair arguments

- 'Ctrl'  
  - Logical-like flag.  
  - If true, the first group X1 is treated as a control group and only comparisons versus this control group are performed (K − 1 comparisons).  
  - If false (default), all pairwise comparisons are performed with Holm-Sidak stepdown adjustment.

- 'Alpha'  
  - Global significance level (scalar in (0,1), default 0.05).

- 'Tail'  
  - Tail of the t-test:
    - 2  → two-sided (default),
    - -1 → left-tailed,
    - 1  → right-tailed.

- 'Display'  
  - Logical-like flag.  
  - If true (default), the function prints a summary of groups and pairwise comparisons to the command window.  
  - If false, the function runs silently and only returns outputs.

## Outputs

The function supports zero, one, or two output arguments.

### No output

holm(X1, X2, ..., XK)  
holm(X1, X2, ..., XK, 'Ctrl', true, 'Alpha', 0.01, 'Display', false)

If no output arguments are requested:

- When Display is true (default), the function prints:
  - a group summary table,
  - degrees of freedom and pooled variance,
  - the pairwise comparison results table (with p-values, Sidak-adjusted alpha, and decision).
- When Display is false, no output is printed and nothing is returned.

### One output

stats = holm(X1, X2, ..., XK, ...)

Returns a structure stats containing group-level summary statistics and global procedure parameters.

### Two outputs

[stats, results] = holm(X1, X2, ..., XK, ...)

Returns both:

- stats   – summary structure
- results – table of pairwise comparisons

### STATS structure

Fields:

- groups – K×1 vector of group indices (1, 2, ..., K), in the order of input arguments.
- N      – K×1 vector of sample sizes per group.
- mean   – K×1 vector of group means.
- std    – K×1 vector of group standard deviations (sample SD).
- df     – degrees of freedom for the pooled variance (Ntot − K).
- s2     – pooled variance across groups.
- alpha  – global significance level used in the procedure.
- tail   – tail of the test (2, -1, or 1).
- k      – number of groups.
- ctrl   – logical flag indicating whether a control group was used.

### RESULTS table

A table with one row per (potential) pairwise comparison. Columns:

- Comparison  – text label indicating the pair of groups (for example "3-1").
- p_value     – p-value of the t-test for that pair, based on the chosen tail.
- Sidak_alpha – Sidak-adjusted alpha used at that step of the Holm-Sidak procedure, or "No comparison made" for skipped comparisons.
- Comment     – textual decision on the null hypothesis, such as:
  - "Reject H0"
  - "Fail to reject H0"
  - "No comparison made" / "H0 is accepted" for comparisons skipped after the first non-significant result.

## Method

All observations from all groups are summarized by group means and standard deviations. Let:

- K be the number of groups,
- N_i the sample size for group i,
- M_i the sample mean for group i,
- S_i the sample standard deviation for group i,
- Ntot = sum_i N_i the total number of observations.

The pooled variance S² is computed as:

S² = Σ_i (N_i − 1) S_i² / (Ntot − K)

with degrees of freedom:

df = Ntot − K

For each pair of groups i and j, the pooled-variance t-statistic is:

t = (M_i − M_j) / sqrt( S² * (1/N_i + 1/N_j) )

The p-value is then derived from the Student's t distribution with df degrees of freedom:

- two-sided:  p = 2 * tcdf(−|t|, df)
- left-tail:  p = tcdf(t, df)
- right-tail: p = tcdf(−t, df)

### Holm-Sidak multiple-comparison control

Let c be the number of comparisons:

- No control group: c = K(K − 1)/2 (all pairwise comparisons).
- With control group: c = K − 1 (comparisons versus group 1 only).

The p-values are sorted in ascending order:

p_(1) ≤ p_(2) ≤ ... ≤ p_(c)

For j = 1, ..., c, the Sidak-adjusted alpha at step j is:

alpha_j = 1 − (1 − alpha)^(1 / (c − j + 1))

The Holm-Sidak stepdown rule is:

1. For j from 1 to c:
   - If p_(j) < alpha_j, reject H0 for the j-th smallest p-value and continue.
   - If p_(j) ≥ alpha_j, fail to reject H0 and stop. All subsequent comparisons are considered "No comparison made" and H0 is accepted for them.

This controls the familywise error rate at the chosen global level alpha.

## Examples

### Example 1 – No control group, printed output

x1 = [7.68 7.69 7.70 7.70 7.72 7.73 7.73 7.76];  
x2 = [7.71 7.73 7.74 7.74 7.78 7.78 7.80 7.81];  
x3 = [7.74 7.75 7.77 7.78 7.80 7.81 7.84];  
x4 = [7.71 7.71 7.74 7.79 7.81 7.85 7.87 7.91];

holm(x1, x2, x3, x4);

The function prints group-level summary statistics, degrees of freedom, pooled variance, and the Holm-Sidak stepdown pairwise comparison table.

### Example 2 – First group as control, alpha = 0.01, silent run with outputs

[stats, results] = holm(x1, x2, x3, x4, 'Ctrl', true, 'Alpha', 0.01, 'Display', false);

No output is printed. stats contains the group summary and procedure parameters, and results contains the pairwise comparison statistics and decisions versus the control.

## Requirements

- MATLAB R14 or later.
- The tcdf function (Statistics and Machine Learning Toolbox in older MATLAB releases; integrated in many newer ones) and basic math functions are required.
- No additional toolboxes are strictly required beyond those providing tcdf.

## Usage Notes

- The legacy syntax using a pooled data vector X and a grouping vector G is not supported in this version. Users are expected to provide one data vector per group.
- Groups are internally labeled as 1, 2, ..., K in the order of the input arguments.
- The function assumes independent groups and pooled variance; it is not intended for repeated-measures or unequal-variance (Welch) designs.

## Citation

If you use this function in a scientific publication, please cite it as:

Cardillo G. (2006). Holm-Sidak t-test: a routine for multiple t-test comparisons. Available on GitHub: https://github.com/dnafinder/holm

## License

The code is provided as-is, without any explicit warranty. Please refer to the repository for licensing details if a LICENSE file is present.
