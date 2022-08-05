function [core_svd, state_pks_full, param] = findSVDensemble(data, param, all_bin_edges)
% Find ensembles using SVD method.
% INPUT:
%     data: N-by-T binary spike matrix, where N is the number of neurons,
%         T is number of frames
%     coords: N-by-2 matrix, containing coordinates of corresponding
%         neurons (for visualization purpose only)
%     param: struct with the following fields:
%         - pks: significant level of spike count per frame, default 4,
%                leave it empty [] if you want an automated threshold
%         - ticut: percentage of cells in a state; default 0.22, leave it
%                empty [] if you want an automated threshold
%         - jcut: another threshold for coactivity: further removes noise;
%                default 0.06
%         - state_cut: maximum number of states allowed, leave it empty if
%                you want automated value
% OUTPUT:
%     core_svd: K-by-1 cell array, where K is the number of identified
%         ensembles, with the indices of core neurons in each ensemble
%     state_pks_full: 1-by-T vector with the active ensemble identity for
%         each frame
%     param: updated parameter structure with the actual values this code
%         used
% 
% For details about this method, see the following paper:
% Carrillo-Reid, et al. "Endogenous sequential cortical activity evoked by
% visual stimuli." Journal of Neuroscience 35.23 (2015): 8813-8828.
% 
% This code is based on Luis Carrillo-Reid's Stoixeon package, but is made
% more efficient and easier to understand.
% Shuting Han, 2018
% 

%% set parameters
if ~isfield(param,'pks'); param.pks = 4; end
if ~isfield(param,'ticut'); param.ticut = 0.22; end
if ~isfield(param,'jcut'); param.jcut = 0.06; end
if ~isfield(param,'state_cut'); param.state_cut = max(floor(size(data, 1) / 2), 1); end
if param.state_cut < 1, error('findSVDensemble:zeroStateCut', 'Cannot continue if no ensembles will be returned'), end

pks = param.pks;
ticut = param.ticut;
jcut = param.jcut;
state_cut = param.state_cut;
if isempty(state_cut); state_cut = round(size(data,1)/4); end

if ~isfield(param, 'findActiveFrames_num_shuff'); param.findActiveFrames_num_shuff = 100; end
findActiveFrames_num_shuff = param.findActiveFrames_num_shuff;
if ~isfield(param, 'findActiveFrames_p'); param.findActiveFrames_p = 0.98; end
findActiveFrames_p = param.findActiveFrames_p;

if ~isfield(param, 'calc_scut_num_shuff'); param.calc_scut_num_shuff = 20; end
calc_scut_num_shuff = param.calc_scut_num_shuff;
if ~isfield(param, 'calc_scut_p'); param.calc_scut_p = 0.98; end
calc_scut_p = param.calc_scut_p;

if ~isfield(param, 'calc_jcut_num_shuff'); param.calc_jcut_num_shuff = 20; end
calc_jcut_num_shuff = param.calc_jcut_num_shuff;
if ~isfield(param, 'calc_jcut_p'); param.calc_jcut_p = 0.99; end
calc_jcut_p = param.calc_jcut_p;

if ~isfield(param, 'SVDStateBinary_fac_cut'); param.SVDStateBinary_fac_cut = 0.4; end
SVDStateBinary_fac_cut = param.SVDStateBinary_fac_cut;
if ~isfield(param, 'SVDStateBinary_p'); param.SVDStateBinary_p = 0.025; end
SVDStateBinary_p = param.SVDStateBinary_p;

if ~isfield(param, 'plot_sequences'); param.plot_sequences = true; end
plot_sequences = param.plot_sequences;
if ~isfield(param, 'plot_ROC'); param.plot_ROC = true; end
plot_ROC = param.plot_ROC;
if ~isfield(param, 'plot_activation_maps'); param.plot_activation_maps = true; end
plot_activation_maps = param.plot_activation_maps;

% Ignore warnings
old_warnings = warning('off', 'stats:pdist:ZeroPoints');


%% find similarity structure
% find high-activity frames
[data_active, pk_indx, pks] = findActiveFrames(data, pks, findActiveFrames_num_shuff, findActiveFrames_p);

% run tf-idf - make this into a function
data_tfidf = calcTFIDF(data_active);
data_tfidf(~isfinite(data_tfidf)) = 0;

% calculate cosine similarity of tf-idf matrix
S_ti = 1 - pdist2(data_tfidf', data_tfidf', 'cosine');

% threshold of noise
if isempty(ticut)
    ticut = calc_scut(data_tfidf, calc_scut_num_shuff, calc_scut_p);
    if ticut >= 1
%         keyboard
        ticut = 0.22; 
        disp(['ticut manually set to ', num2str(ticut)]);
    end % Modified by M.A. to avoid empty similarity after thresholding
end
S_tib = double(S_ti > ticut);
% jaccard similarity, define structures better
if isempty(jcut)
    jcut = calc_jcut(S_tib, calc_jcut_num_shuff, calc_jcut_p);
end
js = 1 - pdist2(S_tib, S_tib, 'jaccard');


%% do SVD, find states and cells
% Find the peaks in the states and the cells in the states
[state_raster, state_pks, SVDStateBinary_fac_cut, max_state_count] = SVDStateBinary(double(js > jcut), state_cut, SVDStateBinary_fac_cut, SVDStateBinary_p);
num_state = size(state_raster,2);

% Set parameters for plotting
n_subplot_columns = ceil(sqrt(num_state));
n_subplot_rows = ceil(num_state/n_subplot_columns);

% get state from full dataset
state_pks_full = zeros(1,size(data,2));
state_pks_full(pk_indx) = state_pks;

if plot_sequences
    % plot
    figure; set(gcf,'color','w')
    raster = zeros(num_state, length(state_pks_full));
    for state = 1:num_state
        raster(state, state_pks_full == state) = 1;
    end
    imagesc(raster)
    colormap(1 - gray)
end
%%  Determine esneble activity in pre and post stim
if num_state > 0
    
    raster = zeros(num_state, length(state_pks_full));
    for state = 1:num_state
        raster(state, state_pks_full == state) = 1;
    end
    activations = raster';
    col_names = cell(1,num_state);
    for i_ens = 1 : num_state
        col_names{i_ens} = ['ensemble_', num2str(i_ens)];
    end
    activation_table = array2table(activations, 'VariableNames', col_names);
    
    all_activations = [activation_table,all_bin_edges];
else
    all_activations = all_bin_edges;

end
max_act_ensemble = find(max(sum(state_raster,1)));

%% find sequences
% find most significant cells for each state
ti_vec = 0.01:0.01:0.1; % cross-validate this threshold
core_svd = cell(num_state,1);
pool_svd = cell(num_state,1);
state_member_raster = zeros(size(data,1), num_state);

if plot_ROC
    figure; clf; set(gcf,'color','w')
    cc = jet(length(ti_vec));
    cc = max(cc-0.3,0);
end

for ii = 1:num_state    
    % pull out all activities in a state
    state_ti_hist = sum(data_tfidf(:, state_pks == ii), 2)';
    state_ti_hist = state_ti_hist / max(state_ti_hist);
    
    % cross-validate ti_cut with cosine similarity
    % the core cells need predict their corresponding state
    if plot_ROC
        subplot(n_subplot_rows,n_subplot_columns,ii); hold on
    end
    auc = zeros(size(ti_vec));
    for n = 1:length(ti_vec)
        core_vec = zeros(size(data_active,1),1);
        core_vec(state_ti_hist > ti_vec(n)) = 1;
        sim_core = 1 - pdist2(data_active', core_vec', 'cosine')';
        [xx, yy, ~, auc(n)] = perfcurve(double(state_pks==ii), sim_core, 1);
        if plot_ROC
            plot(xx,yy,'color',cc(n,:),'linewidth',1);
        end
    end
    if plot_ROC
        plot([0 1],[0 1],'k--')
        xlim([0 1]); ylim([0 1])
        xlabel('FPR'); ylabel('TPR'); title(['ensemble ' num2str(ii)])
    end
    [~, best_indx] = max(auc);
    ti_cut = ti_vec(best_indx);
    state_member_raster(:,ii) = state_ti_hist > ti_cut;
    core_svd{ii} = find(state_member_raster(:, ii));
    pool_svd{ii} = find(state_ti_hist > 0);
end

%% plot core neurons
if plot_activation_maps
    % plot ensemble component cells
    cc_lr = [1 0.8 0.8]; % light red
    cc_r = [1 0.2 0.2]; % dark red
    mksz = 30;
    figure; clf; set(gcf,'color','w')
    for ii = 1:num_state
        subplot(n_subplot_rows,n_subplot_columns,ii); hold on
        scatter(param.coords(:,1),-param.coords(:,2),mksz,'k');
        scatter(param.coords(pool_svd{ii},1),-param.coords(pool_svd{ii},2),mksz,cc_lr,'filled');
        scatter(param.coords(core_svd{ii},1),-param.coords(core_svd{ii},2),mksz,cc_r,'filled');
        title(['ensemble #' num2str(ii)]);
        axis off equal
    end
end

%% update parameters for output
param.pks = pks;
param.ticut = ticut;
param.jcut = jcut;
param.state_cut = state_cut;
param.fac_cut = SVDStateBinary_fac_cut;
param.max_state_count = max_state_count;
param.max_act_ensemble = max_act_ensemble;
param.all_activations = all_activations;

disp(['found ', num2str(length(core_svd)), ' ensembles'])

% Restore warning state
warning(old_warnings)
