function [stats, results] = holm(varargin)
%HOLM Holm-Sidak procedure for multiple Student's t-tests.
%
%   HOLM(X1, X2, ..., XK)
%   HOLM(X1, X2, ..., XK, 'Ctrl', CTRL, 'Alpha', ALPHA, 'Tail', TAIL, 'Display', DISPLAY)
%
%   Performs the Holm-Sidak stepdown procedure for multiple comparisons
%   between K independent groups using pooled-variance Student's t-tests.
%
%   Each Xi is the data vector for the i-th group (row or column, with
%   possibly different lengths).
%
%   Name–Value pair arguments:
%     'Ctrl'    - logical-like flag (true/false, 0/1, 'on'/'off',
%                 'true'/'false', 'yes'/'no')
%                 true  -> the first group (X1) is a control group
%                 false -> no control group (default)
%
%     'Alpha'   - significance level (scalar in (0,1), default = 0.05)
%
%     'Tail'    - tail of the test:
%                 2  -> two-sided (default)
%                -1  -> left-tailed
%                 1  -> right-tailed
%
%     'Display' - logical-like flag controlling command-window output:
%                 true  -> print tables and messages (default)
%                 false -> run silently, only return outputs
%
%   Outputs:
%     stats   - structure with summary statistics for groups and procedure
%     results - table with pairwise comparison results
%
%   If no output arguments are requested, the function behaves in a
%   MATLAB-style way: it prints to the command window (unless Display is
%   set to false) and does not return stats or results.
%
%   Example:
%
%   x1 = [7.68 7.69 7.70 7.70 7.72 7.73 7.73 7.76];
%   x2 = [7.71 7.73 7.74 7.74 7.78 7.78 7.80 7.81];
%   x3 = [7.74 7.75 7.77 7.78 7.80 7.81 7.84];
%   x4 = [7.71 7.71 7.74 7.79 7.81 7.85 7.87 7.91];
%
%   % No control group, with printed output:
%   holm(x1, x2, x3, x4)
%
%   % First sample is the control group, alpha = 0.01, no printed output:
%   [stats, results] = holm(x1, x2, x3, x4, 'Ctrl', true, 'Alpha', 0.01, 'Display', false);
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2006). Holm-Sidak t-test: a routine for multiple t-test
%   comparisons. Available on GitHub: https://github.com/dnafinder/holm

% -------------------------------------------------------------------------
% Input parsing: data arguments and Name–Value pairs
% -------------------------------------------------------------------------

narginchk(2, Inf);

% Separate data arguments (groups) from Name–Value pairs
isNameArg = @(c) ischar(c) || (isstring(c) && isscalar(c));
firstNV   = find(cellfun(isNameArg, varargin), 1, 'first');

if isempty(firstNV)
    dataArgs = varargin;
    nvArgs   = {};
else
    dataArgs = varargin(1:firstNV-1);
    nvArgs   = varargin(firstNV:end);
end

k = numel(dataArgs);
assert(k >= 2, 'HOLM requires at least two groups.');

% Parse Name–Value pairs
p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Ctrl',    false, @validateLogicalLike);
addParameter(p, 'Alpha',   0.05,  @(x) validateattributes(x,{'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}, mfilename, 'Alpha'));
addParameter(p, 'Tail',    2,     @(x) isnumeric(x) && isscalar(x) && ismember(x,[-1 1 2]));
addParameter(p, 'Display', true,  @validateLogicalLike);
parse(p, nvArgs{:});

ctrl        = logical(normalizeLogicalLike(p.Results.Ctrl));
alpha       = p.Results.Alpha;
tail        = p.Results.Tail;
displayFlag = logical(normalizeLogicalLike(p.Results.Display));

% Collect data: N, mean, std for each group
N  = zeros(k,1);
Md = zeros(k,1);
Sd = zeros(k,1);

for i = 1:k
    Xi = dataArgs{i};
    validateattributes(Xi, {'numeric'}, ...
        {'vector','real','finite','nonnan','nonempty'}, ...
        mfilename, sprintf('X%d',i), i);
    Xi   = Xi(:);          % column vector
    N(i) = numel(Xi);
    Md(i)= mean(Xi);
    Sd(i)= std(Xi);        % sample standard deviation
end

groups = (1:k).';

% -------------------------------------------------------------------------
% Main computation
% -------------------------------------------------------------------------

n  = sum(N);
df = n - k;                            % degrees of freedom
s2 = sum((N-1).*Sd.^2) / df;           % pooled variance

% Group summary table
Tgroups = table(groups, N, Md, Sd, ...
    'VariableNames', {'Group','N','Mean','Standard_deviation'});

% Determine number of comparisons
if ~ctrl
    % No control group: all pairwise comparisons
    c = k*(k-1)/2;
else
    % With control group: only comparisons vs group 1
    c = k-1;
end

pVals = nan(c,1);
pb    = cell(c,4); % {Comparison, p_value, Sidak_alpha/No comparison, Comment}
count = 0;

% Compute t and p for each comparison
if ~ctrl
    % All pairwise comparisons
    for I = 1:(k-1)
        for J = (I+1):k
            count = count + 1;
            t = (Md(I) - Md(J)) / sqrt(s2 * sum(1./N([I J])));
            switch tail
                case -1 % left-tailed
                    pVals(count) = tcdf(t, df);
                case 1  % right-tailed
                    pVals(count) = tcdf(-t, df);
                case 2  % two-sided
                    pVals(count) = 2*tcdf(-abs(t), df);
            end
            pb(count,1:2) = {sprintf('%d-%d', I, J), pVals(count)};
        end
    end
else
    % Comparisons vs control (group 1)
    I = 1;
    for J = 2:k
        count = count + 1;
        t = (Md(I) - Md(J)) / sqrt(s2 * sum(1./N([I J])));
        switch tail
            case -1 % left-tailed
                pVals(count) = tcdf(t, df);
            case 1  % right-tailed
                pVals(count) = tcdf(-t, df);
            case 2  % two-sided
                pVals(count) = 2*tcdf(-abs(t), df);
        end
        pb(count,1:2) = {sprintf('%d-%d', I, J), pVals(count)};
    end
end

% Sort p-values for Holm-Sidak stepdown
[pSorted, idx] = sort(pVals);
pb             = pb(idx,:);

% Holm-Sidak adjusted alpha values
cIdx      = (1:c).';
alphaEff  = 1 - (1 - alpha).^(1 ./ (c - cIdx + 1));

% Stepdown decision
comp = true; % comparison checker
for j = 1:c
    if comp
        if pSorted(j) < alphaEff(j)
            pb(j,3:4) = {alphaEff(j), 'Reject H0'};
        else
            pb(j,3:4) = {alphaEff(j), 'Fail to reject H0'};
            comp = false; % no further comparison is required
        end
    else
        pb(j,3:4) = {'No comparison made', 'H0 is accepted'};
    end
end

% Pairwise comparison results table
results = cell2table(pb, 'VariableNames', {'Comparison','p_value','Sidak_alpha','Comment'});

% -------------------------------------------------------------------------
% Build outputs
% -------------------------------------------------------------------------

stats = struct();
stats.groups      = groups;
stats.N           = N;
stats.mean        = Md;
stats.std         = Sd;
stats.df          = df;
stats.s2          = s2;
stats.alpha       = alpha;
stats.tail        = tail;
stats.k           = k;
stats.ctrl        = ctrl;

% Display (if requested)
if displayFlag
    disp('HOLM-SIDAK PROCEDURE FOR MULTIPLE STUDENT''S T TESTS');
    disp(' ');
    disp(Tgroups);
    fprintf('Degrees of freedom: %d - Combined variance: %0.4f\n', df, s2);
    disp(' ');
    disp(results);
end

% If no output is requested, clear outputs (MATLAB-style)
if nargout == 0
    clear stats results
end

end

% -------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function tf = validateLogicalLike(x)
%VALIDATELOGICALLIKE Helper for inputParser: check logical-like values.
    try
        normalizeLogicalLike(x);
        tf = true;
    catch
        tf = false;
    end
end

function y = normalizeLogicalLike(x)
%NORMALIZELOGICALLIKE Convert various logical-like inputs to true/false.
    if islogical(x)
        y = x;
    elseif isnumeric(x) && isscalar(x)
        y = (x ~= 0);
    elseif ischar(x) || (isstring(x) && isscalar(x))
        s = lower(char(x));
        if any(strcmp(s, {'true','on','yes'}))
            y = true;
        elseif any(strcmp(s, {'false','off','no'}))
            y = false;
        else
            error('Invalid logical-like value: %s', s);
        end
    else
        error('Invalid type for logical-like option.');
    end
end
