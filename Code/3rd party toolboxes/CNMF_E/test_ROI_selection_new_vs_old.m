% compare outputs of CNMF-e and old protocol
spikes_filename = get_filename_of('spikes', animal_ID);
spikes = load_variable(spikes_filename, 'spikes');
dff_filename = get_filename_of('dFF', animal_ID);
dFF = load_variable(dff_filename, 'dFF');
%% Get ROI info
ROI_info_filename = get_filename_of('ROI_info', animal_ID);
% load ROis
ROIs_old = load_variable(ROI_info_filename, 'ROI_info');
% position 5$
locs = NaN(length(ROIs_old),2);
for ir = 1:length(ROIs_old)
    iloc = nanmean(cell2mat(ROIs_old(ir,5)));
    locs (ir,:) = iloc;
end

%% Plot ROIs old
fig_fov = figure('pos', [10 10 1600 800]);
subplot(1,2,1)
title('Manual segmentation')
imshow(mean_images) % this comes form CNMF-e main script
hold on
for ir = 1:length(ROIs_old)
    pos = [locs(ir,1), locs(ir,2)];
    scatter(pos(1), pos(2))
    text(locs(ir,1), locs(ir,2),num2str(ir));
end
hold off

% plot new
subplot(1,2,2)
title('CNMF-e')
Cn = (mean_images); % this comes form CNMF-e main script
Coor = neuron_concat.Coor;
Aor = neuron_concat.A;
CC = plot_contours(Aor, Cn,[], 1, [], Coor);

%% get example cell
old_cell = 1;
new_cell = 14;

% get the spikes & dFF
sp_old = spikes;
all_spikes_old = cell2mat(sp_old(:,:));
spikes_old = all_spikes_old(old_cell,:);
dFF_old = dFF(old_cell,:);
dFF_old = cell2mat(dFF_old);
all_dff_old = cell2mat(dFF(:,:));


% Get the CNFM-e spikes & dFF
all_spikes_new =  neuron_concat.S;
spikes_new = neuron_concat.S(new_cell,:);
all_dFF_new = neuron_concat.C_raw;
dFF_new = all_dFF_new(new_cell,:);

% Plot spikes
fig_spikes= figure('pos', [10 10 1600 800]);
plot(spikes_old/max(spikes_old), 'color', 'b' ,'LineWidth', 2)
hold on
plot(spikes_new/max(spikes_new), 'color', 'r','LineWidth', 2)
legend({'old', 'new'})
title('SPIKES')

% Plot dFF
% fig_dFF= figure('pos', [10 10 1600 800]);
plot(dFF_old/max(dFF_old),'color', 'b', 'LineWidth', .3, 'LineStyle','--')
hold on
plot(dFF_new/max(dFF_new),'color', 'r','LineWidth', .3, 'LineStyle','--')
legend({'old', 'new'})
title('Fluorescence')

%% calculate correlations
% generate a matrix to find highly correlated cells between both methods
n_old = size(all_spikes_old,1);
n_new = size(all_spikes_new,1);

% DFF:  generate a matrix
XC = NaN(n_old,n_new);
for i_old = 1:n_old
    this_data_old = all_dff_old(i_old, :)/max(all_dff_old(i_old, :));
    for i_new = 1:n_new
        % allocate
         this_data_new = all_dFF_new(i_new, :)/ max(all_dFF_new(i_new, :));
         this_xr = corrcoef(this_data_old, this_data_new);
         XC(i_old, i_new) = this_xr(1,2);
    end
end
%%
figure, 
imagesc(XC>0.3)
caxis([0 0.25])
ylabel('Old cells')
xlabel('New cells')




