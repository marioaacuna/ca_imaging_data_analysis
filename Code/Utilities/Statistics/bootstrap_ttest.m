function [p, t, df] = bootstrap_ttest(DATA, LABELS, pairedTest, permutationsNo, withReplacement, showFigure)
% TTEST.BOOTSTRAP performs a randomization test through bootstrapping.

if ~exist('permutationsNo','var') || isempty(permutationsNo), permutationsNo = 10^4; end
if ~exist('withReplacement','var') || isempty(withReplacement), withReplacement = true; end
if ~exist('showFigure','var') || isempty(showFigure), showFigure = false; end

if ~pairedTest  % linearize inputs for independent samples t-test
    DATA = DATA(:);
    LABELS = LABELS(:);
    [~, ~, LABELS] = unique(LABELS, 'stable');
else
    LABELS = ones(size(DATA,1),1);
end
% Remove NaNs
rows = any(isnan(DATA'),1)';
DATA(rows,:) = [];
LABELS(rows,:) = []; 

% Generate a new randomization seed; be aware that following line makes each run unique
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% --- Two-sample paired t-test ---
if pairedTest
    % This test assumes a pre-post design with a subject per row.
    % If we shuffle subjects, we lose information about within-subject variance,
    % hence what we want to randomize is actually the pre/post label.
    % For computational ease, the sign of the pre-post difference will be randomly
    % assigned to either 1 or -1.
    
    % Estimate differences between pre- and post-treatment values
    x_pre  = DATA(:,1);
    x_post = DATA(:,2);
    x = x_pre - x_post;
    LABELS = sign(x);  % These are the values that will be shuffled
    n = length(x);
    Estimate0 = mean(x); % This is our sample case. H1: a pre-post effect exists and it has this mean value
    
    if withReplacement
        % Create random [-1 +1] values
        BS_Labels = randi(2,n,permutationsNo);
        BS_Labels = (BS_Labels -1.5) .*2;
    else
        % Shuffle labels, keeping constant the number of +1s and -1s
        BS_Labels = NaN(n,permutationsNo);
        parfor p = 1:permutationsNo
            BS_Labels(:,p) = LABELS(randperm(n));
        end
    end
    
    % --- do bootstrap ---
    abs_x = abs(x); % Remove +/- sign
    X = bsxfun (@times, BS_Labels, abs_x); % Multiply +1/-1s to the absolute differences
    Estimate_BS = mean(X); % bootstrapped estimation of differences. H0: there is no effect, the mean is 0

    % Calculate degrees of freedom and t-value like in a parametric t-test
%     df = n -1;
%     t = mean(x) / (std(x)/sqrt(n));

else
    % --- Two-sample un-paired t-test ---
    % With this test, 2 independent samples need to be tested. The null-hypothesis
    % is that there is no difference between sample means. Labels of all samples will be shuffled.
    % For ease of comparison with a parametric t-test, the Student-s t distribution will be 
    % used to determine statistical significance.
    
    % Estimate t-score from the sample
    x0 = DATA(LABELS==1);
    x1 = DATA(LABELS==2);
    n0 = length(x0);
    n1 = length(x1);
    Estimate0 = (mean(x0)-mean(x1));% / sqrt((var(x0)/n0) + (var(x1)/n1));
    
    if withReplacement
        BS_Labels = randi(2,n0+n1,permutationsNo);
    else
        BS_Labels = NaN(n0+n1,permutationsNo);
        parfor p = 1:permutationsNo
            BS_Labels(:,p) = LABELS(randperm(n0+n1));
        end
    end
    
    % Set H0 as true (i.e. mean is 0)
    GM = mean(DATA); % grand mean
    DATA2 = [x0-mean(x0); x1-mean(x1)] +GM;

    % --- do bootstrap ---
    Estimate_BS = NaN(permutationsNo,1);
    parfor loop = 1:permutationsNo
        thisLabels = BS_Labels(:,loop);
        x0 = DATA2(thisLabels==1);
        x1 = DATA2(thisLabels==2);
        
        Estimate_BS(loop) = (mean(x0)-mean(x1));% / sqrt((var(x0)/length(x0)) + (var(x1)/length(x1)));
    end
    
    % Calculate degrees of freedom with Satterthwaite's approximation (from ttest2, with unequal variances)
%     t = Estimate0;
%     s2x = var(x0); s2y = var(x1);
%     s2xbar = s2x ./n0;  s2ybar = s2y ./n1;
%     df = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (n0-1) + s2ybar.^2 ./ (n1-1));
end
t=NaN; df=NaN;

% Compare original estimate with bootstraped distribution. Estimate a 2-tailed p-value
p = (sum(abs(Estimate_BS) >= abs(Estimate0)) +1) / (permutationsNo+1);


% Show a figure, if required
if showFigure
    figure('color','w');
    clf; hold on; Lh=[];
    [H,x] = hist(Estimate_BS,100);
    H = H ./permutationsNo .*100;
    Bh = bar(x,H,'hist');
    set(Bh,'Facecolor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);
    m = mean(Estimate_BS);
    Lh(2) = line([m m],[0 permutationsNo],'color','k','linewidth',2,'YLimInclude','off');
    Lh(1) = line([Estimate0 Estimate0],[0 permutationsNo],'color','r','linewidth',3,'YLimInclude','off');
    uistack(Lh,'top')

    legend(Lh,{'sample','bootstraped mean'});
    xlabel('bootstraped distribution'); ylabel('bootstrap samples [%]');
    title(['\itp\rm: ' num2str(p) ' (' num2str(permutationsNo,'%1.1g') ' permutations)']);
    set(gca,'Layer','top','tickdir','out');
end 
