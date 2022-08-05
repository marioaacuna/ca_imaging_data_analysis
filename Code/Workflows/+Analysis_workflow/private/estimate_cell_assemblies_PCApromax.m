function RESULTS = estimate_cell_assemblies_PCApromax(data, test_swap_cells, verbose)
%ESTIMATE_CELL_ASSEMBLIES_PCAPROMAX Summary of this function goes here
%   Detailed explanation goes here

if ~exist('verbose','var'), verbose = true; end
if ~exist('test_swap_cells','var'), test_swap_cells = false; end

% Unpack user inputs
if iscell(data)
    data = cell2mat(data);
end

% Read general_configs and LOGGER
global GC LOGGER
correlation_type = GC.assemblies_correlation_type;
correlation_merge_assemblies = GC.correlation_merge_assemblies;
n_permutations_significance = GC.assemblies_n_permutations_significance;
zmax_cutoff = GC.assemblies_zmax_cutoff;

% Initialize output variable
RESULTS = struct('z_max',NaN, 'n_silent_cells',NaN, 'assemblies',[], 'label_permutation_p',NaN, ...
    'cell_swap_p',NaN, 'strong_assemblies',[], 'weak_assemblies',[], 'cells_not_in_assemblies',[]);

% Remove cells that never fired
n_cells = size(data, 1);
silent_cells = find(sum(data, 2) == 0);
data(silent_cells, :) = [];
RESULTS.n_silent_cells = length(silent_cells);

% Get number of active cells and frames
[n_active_cells, n_frames] = size(data);

% Transform the response to z-score    
Q = zscore(data, [], 2);
% Perform PCA
[PCs, ~, eigenvals] = pca(Q', 'Algorithm','svd', 'Economy',false);

% Dataset dimensionality reduction is obtained by only keeping principal
% components (PCs) with eigenvalues greater than lambda_max, a theoretical
% lower bound to the eigenvalues of informative PCs given by the Marcenko
% Pastur distribution:
% lambda = (1 + sqrt(N / T) ^2) + N^(-2/3)
lambda = (1 + sqrt(n_active_cells / n_frames)) ^ 2;
% Fix this value by adding the Tracy-Widom factor
tracy_widom_correction = n_active_cells ^ (-2/3);
lambda_max = lambda + tracy_widom_correction;
% The PCs larger that lambda_max are termed "signal" components
signal_PCs = find(eigenvals > lambda_max);
n_signal_PCs = length(signal_PCs);
% If no signal PCs survived, quit algorithm
if n_signal_PCs < 1
    return
end

% Perform promax PC rotation
if signal_PCs == 1
    rPCs = PCs(:, signal_PCs);
else
    try
        [rPCs, ~] = rotatefactors(PCs(:, signal_PCs), 'Method','promax', 'Maxit',100000);
    catch ME
        if strcmp(ME.identifier, 'stats:rotatefactors:IterationLimit')
            if verbose
                LOGGER.error('Promax did not converge. Impossible to determine assemblies. Skipping dataset')
            end
            return
        else
            rethrow(ME)
        end
    end
end

% Normalize and flip axes such that largest value on any PCs has a positive value
for ipc = 1:n_signal_PCs
    rPCs(:,ipc) = rPCs(:,ipc) / norm(rPCs(:,ipc));
    if max(rPCs(:,ipc)) <= abs(min(rPCs(:,ipc)))
        rPCs(:,ipc) = -rPCs(:,ipc);
    end
end
% Compute projections of each cell onto the rotated PCs
z_rPCs = zscore(rPCs);
rPC_proj = max(z_rPCs, [], 2);
z_max = fit_zmax_to_data(rPC_proj, zmax_cutoff);
RESULTS.z_max = z_max;

% Reconstruct putative cell assemblies
cell_assemblies          = {};
cell_assemblies_vector   = zeros(size(rPCs));
cell_assemblies_strength = zeros(n_signal_PCs, 1);
PCs_to_del = [];
for ipc = 1:n_signal_PCs
    % Sort values
    [values, inds] = sort(z_rPCs(:,ipc), 'descend');
    % Find out which values are above z_max
    above_z_max = values >= z_max;
    if ~any(above_z_max)
        PCs_to_del = [PCs_to_del; ipc];
    else
        % Get indices of cells in thie assembly
        cells_in_assembly = inds(above_z_max)';
        cell_assemblies{end+1} = cells_in_assembly;
        % Compute the strength of the assembly
        cell_assemblies_vector(cells_in_assembly, ipc) = rPCs(cells_in_assembly, ipc);
        cell_assemblies_strength(ipc) = norm(cell_assemblies_vector(:,ipc));
        cell_assemblies_vector(:,ipc) = cell_assemblies_vector(:,ipc) / cell_assemblies_strength(ipc);
    end
end
% Remove PCs without assemblies on them
if ~isempty(PCs_to_del)
    rPCs(:,PCs_to_del)                   = [];
    cell_assemblies_vector(:,PCs_to_del) = [];
end

% Check for assemblies that are very similar and merge them
rotPCsMerged = rPCs;
cell_assemblies_vectorMerged = cell_assemblies_vector;
cell_assembliesMerged = cell_assemblies;
n_iter = 0;
% We calculate the projection between the assemblyVectors and if they are
% bigger than similarityThresh we merge the assembly pair
while true
    % Projection
    similarityAssemblies = cell_assemblies_vectorMerged' * cell_assemblies_vectorMerged;
    similarityAssemblies(logical(triu(similarityAssemblies))) = 0;
    % Check pairs too similars, and delete them one at a time
    [similar1, similar2] = find(similarityAssemblies >= correlation_merge_assemblies);
    if isempty(similar1)
        break
    end
    vals = zeros(size(similar1));
    for isim = 1:length(similar1)
        vals(isim) = similarityAssemblies(similar1(isim), similar2(isim));
    end
    [~, indMax] = max(vals);
    % the merged assemblies is the (thresholded) vectorial sum of their vectors
    summedVectors = (rotPCsMerged(:,similar1(indMax)) + rotPCsMerged(:,similar2(indMax)));
    rotPCsMerged(:,[similar1(indMax), similar2(indMax)]) = [];
    cell_assemblies_vectorMerged(:,[similar1(indMax), similar2(indMax)]) = [];
    cell_assembliesMerged([similar1(indMax), similar2(indMax)]) = [];
    newVect = summedVectors / norm(summedVectors);
    rotPCsMerged(:,end+1) = newVect;
    cell_assembliesMerged{end+1} = find(zscore(newVect) >= z_max)';
    newVect(zscore(newVect) < z_max) = 0;
    cell_assemblies_vectorMerged(:,end+1) = newVect / norm(newVect);
    n_iter = n_iter + 1;
end
% Re-copy values in older variables
rPCs = rotPCsMerged;
cell_assemblies = cell_assembliesMerged;
clear cell_assemblies_vectorMerged rotPCsMerged cell_assembliesMerged;
if n_iter > 0
    if verbose
        LOGGER.trace([num2str(n_iter), ' assemblies were merged because correlation >= ', num2str(correlation_merge_assemblies)])
    end
end

% Remove assemblies that are left with only one cell
n_cells_in_assembly = cellfun(@length, cell_assemblies);
to_delete = find(n_cells_in_assembly < 2);
cell_assemblies(to_delete) = [];
rPCs(:, to_delete)=[];
n_signal_PCs = size(rPCs,2);
if ~isempty(to_delete)
    if verbose
        LOGGER.trace([num2str(length(to_delete)), ' assemblies were removed because they had only one cell left in them'])
    end
end
if n_signal_PCs == 0
    return
end

% Compute the pairwise cell activity correlation matrix
C = corr(Q', 'type',correlation_type);

% Assess assemblies significance by calculating correlation
correlation_in_assemblies = NaN(n_signal_PCs, 1);
p_random_cells = NaN(n_signal_PCs, 1);
p_swap_cells = cell(n_signal_PCs, 1);

for i_assembly = 1:length(cell_assemblies)
    n_cells_in_assembly = length(cell_assemblies{i_assembly});
    base_idx = nchoosek(1:n_cells_in_assembly, 2);
    n_base_idx = size(base_idx, 1);

    % Compute correlation in the assembly
    idx = nchoosek(sort(cell_assemblies{i_assembly}), 2);
    idx_linear = reshape(sub2ind(size(C), idx(:,1), idx(:,2)), n_base_idx, []);
    correlation_in_assemblies(i_assembly) = mean(C(idx_linear));

    % Compute null distribution of correlation coefficients for a random
    % assembly with the same number of cells
    null_distribution = NaN(n_permutations_significance, 1);
    parfor irep = 1:n_permutations_significance
        shuffledAssembly = randi(n_active_cells, 1 ,n_cells_in_assembly);
        idx = nchoosek(sort(shuffledAssembly), 2);
        idx_linear = reshape(sub2ind(size(C), idx(:,1), idx(:,2)), n_base_idx, []);
        null_distribution(irep) = mean(C(idx_linear));
    end
    % Compute p-value
    p_random_cells(i_assembly) = (sum(null_distribution(:) >= correlation_in_assemblies(i_assembly)) +1) / (n_permutations_significance+1);
    
    % Compute null distribution of correlation coefficients within the assembly
    % when swapping up to n-1 cells, where n is the number of cells in the
    % assembly
    if test_swap_cells
        cells_not_in_this_assembly = setdiff(1:n_active_cells, cell_assemblies{i_assembly});
        cells_not_in_this_assembly = cells_not_in_this_assembly(:);
        n_cells_not_in_this_assembly = length(cells_not_in_this_assembly);
        p = NaN(n_cells_in_assembly-1, 2);
        ii = 1;
        for n_cells_to_keep = n_cells_in_assembly-1:-1:1
            cells_in_surrogate_assembly = nchoosek(1:n_cells_in_assembly, n_cells_to_keep);
            n_combinations = size(cells_in_surrogate_assembly, 1);
            p_residual_correlation = NaN(n_combinations, 1);
            parfor i_comb = 1:n_combinations
                % Add random cells to complete the assembly
                cells_from_original_assembly = cell_assemblies{i_assembly}(cells_in_surrogate_assembly(i_comb, :));
                n_cells_to_add = n_cells_in_assembly - n_cells_to_keep;
                new_assembly = [repmat(cells_from_original_assembly, n_permutations_significance, 1), cells_not_in_this_assembly(randi(n_cells_not_in_this_assembly, n_permutations_significance, n_cells_to_add))];
                base_idx = nchoosek(1:n_cells_in_assembly, 2);
                n_base_idx = size(base_idx, 1);
                idx = NaN(n_base_idx * n_permutations_significance, 2);
                for irep = 1:n_permutations_significance
                    idx((irep-1)*n_base_idx+1:(irep-1)*n_base_idx+n_base_idx, :) = reshape(new_assembly(irep, base_idx), [], 2);
                end
                idx_linear = reshape(sub2ind(size(C), idx(:,1), idx(:,2)), n_base_idx, []);
                null_distribution = C(idx_linear);
                % Compute average p-value
                p_residual_correlation(i_comb) = mean((sum(null_distribution >= correlation_in_assemblies(i_assembly), 1) +1) / (n_permutations_significance+1));
            end
            p(ii, 1) = n_cells_to_keep;
            p(ii, 2) = mean(p_residual_correlation);
            ii = ii + 1;
        end
        p_swap_cells{i_assembly} = array2table(p, 'VariableNames',{'n_cells_kept','p'});
    end
end

% Make table with cell and assembly ID
all_cell_idx = (1:n_cells)';
all_cell_idx(silent_cells) = [];
active_cells = 1:n_active_cells;
cells_not_in_assemblies = setdiff(active_cells, cell2mat(cell_assemblies));

RESULTS.assemblies = [];
RESULTS.label_permutation_p = p_random_cells;
RESULTS.cell_swap_p = p_swap_cells;
RESULTS.strong_assemblies = find(p_random_cells <= 0.05);
RESULTS.weak_assemblies = find(p_random_cells > 0.05);
RESULTS.cells_not_in_assemblies = cells_not_in_assemblies;

for i_assembly = 1:length(cell_assemblies)
    cells_in_assembly = cell_assemblies{i_assembly};
    % Compute mean correlation of each cell to rest of the assembly
    mean_correlation_within_assembly = NaN(length(cells_in_assembly), 1);
    for i_cell = 1:length(cells_in_assembly)
        idx = cells_in_assembly;
        idx(i_cell) = [];
        mean_correlation_within_assembly(i_cell) = mean(C(cells_in_assembly(i_cell), idx));
    end
    % Compute mean correlation of each cell to cells in other assemblies
    mean_correlation_other_assemblies = NaN(length(cells_in_assembly), 1);
    for i_cell = 1:length(cells_in_assembly)
        for i_assembly_inner = 1:length(cell_assemblies)
            if i_assembly_inner == i_assembly
                continue
            end
            cells_in_other_assembly = cell_assemblies{i_assembly_inner};
            mean_correlation_other_assemblies(i_cell) = mean(C(cells_in_assembly(i_cell), cells_in_other_assembly));
        end
    end
    % Compute mean correlation of each cell to cells that did not fall in any
    % other assembly
    mean_correlation_outside_assemblies = NaN(length(cells_in_assembly), 1);
    for i_cell = 1:length(cells_in_assembly)
        mean_correlation_outside_assemblies(i_cell) = mean(C(cells_in_assembly(i_cell), cells_not_in_assemblies));
    end

    % Store values
    true_cell_idx = all_cell_idx(cells_in_assembly);
    res = [true_cell_idx(:), repmat(i_assembly, length(cells_in_assembly), 1), mean_correlation_within_assembly(:), mean_correlation_other_assemblies(:), mean_correlation_outside_assemblies(:)];
    RESULTS.assemblies = [RESULTS.assemblies; res];
end
RESULTS.assemblies = array2table(RESULTS.assemblies, 'VariableNames',{'cell', 'assembly', 'mean_correlation_within_assembly', 'mean_correlation_other_assemblies', 'mean_correlation_outside_assemblies'});

if verbose
    LOGGER.trace(['Identified ', num2str(length(RESULTS.strong_assemblies)), ' cell assemblies'])
end

%% MLint exceptions
%#ok<*AGROW,*PFBNS>
