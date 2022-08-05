function [H, corrP] = adjust_Pvalues(p, alpha0, method)
% ADJUST_PVALUES perform correction of p-values due to multiple comparisons bias, using
% the step-wise Holm procedure. It supports both 'Sidak' and 'Bonferroni' algorhitms.
% 
% Further details:
% http://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_holms_multiple_comparison_test.htm
% https://scholar.google.ch/citations?view_op=view_citation&hl=en&user=Tyd1NtYAAAAJ&citation_for_view=Tyd1NtYAAAAJ:DIubQTN3OvUC
% http://pages.infinit.net/rlevesqu/Syntax/Unclassified/AdjustedP-ValuesAlgorithms.txt
% 

if ~exist('alpha0','var') || isempty(alpha0), alpha0 = 0.05; end
if ~exist('method','var') || isempty(method), method = 'Sidak'; end

% Convert to array
if iscell(p)
    input_is_cell = true;
    p = cell2mat(p);
else
    input_is_cell = false;
end

%% Apply correction
P = [(1:length(p))' p(:)];
P2 = sortrows(P, 2); % ascending order
IDX = ~isfinite(P2(:,2));
P = P2(~IDX,:);
C = size(P,1);
switch lower(method(1))
    case 's' % Sidak
        P(:,3) = 1 - (1-P(:,2)).^(C-(1:C)'+1);
    case 'b' % Bonferroni
        P(:,3) = P(:,2) .*C;
end

%% Holm procedure for selection
% Check comparisons
P(:,4) = P(:,3) <= alpha0;
% Find first accepted H0
acc = find(P(:,4)==0, 1, 'first');
if ~isempty(acc)
    P(acc:end,4) = 0; % all subsequent H0 hypotheses must be accepted
end

% Generate output
P = sortrows(P, 1); % Re-sort according to original order
% Add back non-finite p-values
p = NaN(length(IDX), size(P, 2));
p(P(:, 1), :) = P;
P = p;
P(:,1) = 1:size(P,1);
P(isnan(P(:, 2)), 2:4) = repmat([NaN, 1, 0], sum(IDX), 1);

H = logical(P(:,4)); % Do output a logical vector that answers the question "Does this p-value cross the significance threshold?"
corrP = P(:,3);
corrP(corrP > 1) = 1;

% Reconvert ouput
if input_is_cell
    corrP = num2cell(corrP);
end
