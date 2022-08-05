%           |  ALL   |  SIG   |
%           | LP SAL | LP SAL |
% __________|________|________|
% SP  |CCI  | 599 56 | 181 29 |
% ____|sham | 353 30 | 112  8 |
% HPS |CCI  | 585 56 | 181 29 |
%     |sham | 340 29 | 112  8 |




x = linspace(0, 600, 5100);




sig = cell(2, 2);
for i = 1:2
    for ii = 1:2
        a = results.HPS.DATA_p{i, ii}(:, 6) + results.HPS.DATA_p{i, ii}(:, 7);
        b = a > 0;
        sig{i, ii} = [results.HPS.DATA_p{i, ii}(:,1:2) b];
    end
end

sig = cell(2, 2);
sig_all = cell(2, 2);
for i = 1:2
    for ii = 1:2
        sig{i, ii} = sig_all{i, ii}(any(sig_all{i, ii}(:, 3), 3),:);
        sig{i, ii}(:, end) = [];
    end
end




sig_all = cell(2, 2);
for i = 1:2
    for ii = 1:2
        a = results.HPS.DATA_p{i, ii}(:, 6) + results.HPS.DATA_p{i, ii}(:, 7);
        b = a > 0;
        sig_all{i, ii} = [results.HPS.DATA_p{i, ii}(:,1:2) b];
    end
end





sig = cell(2, 2)
for i = 1:2
    for ii = 1:2
        sig{i, ii} = sig_all{i, ii}(any(sig_all{i, ii}(:, 3), 3),:);
        sig{i, ii}(:, end) = [];
    end
end



sig = cell(2, 2)
for i = 1:2
    for ii = 1:2
        sig{i, ii} = sig_all{i, ii}(any(sig_all{i, ii}(:, 3), 3),:);
        sig{i, ii}(:, 3) = ones(numel(sig{i, ii}(:, 1)), 1);
    end
end
sig2 = [sig(1, :)'; sig(2, :)'];



sig = cell(2, 2)
for i = 1:2
    for ii = 1:2
        sig{i, ii} = sig_all{i, ii}(any(sig_all{i, ii}(:, 3), 3),:);
        if i == 1 && ii == 1
            sig{i, ii}(:, 3) = ones(numel(sig{i, ii}(:, 1)), 1);
        elseif i == 1 && ii == 2
            sig{i, ii}(:, 3) = 2.*ones(numel(sig{i, ii}(:, 1)), 1);
        elseif i == 2 && ii == 1
            sig{i, ii}(:, 3) = 3.*ones(numel(sig{i, ii}(:, 1)), 1);
        elseif i == 2 && ii == 2
            sig{i, ii}(:, 3) = 4.*ones(numel(sig{i, ii}(:, 1)), 1);
        end
    end
end
sig2 = [sig(1, :)'; sig(2, :)'];
sig3 = {cat(1, sig2{:})}


