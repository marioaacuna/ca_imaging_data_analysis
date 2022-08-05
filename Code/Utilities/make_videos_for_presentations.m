% Create a video of FOV and traces from different sessions for presentations

%% Clear variables
clear
close all
clc
global GC

% %% Load CNMFe files
% 
% cnmfe_folder = 'C:\Users\acuna\Documents\CNMF_E';
% cnmfe_loaded = 0;
% if ~exist('cnmfe_loaded', 'var') || ~cnmfe_loaded
%     addpath(fullfile(cnmfe_folder, 'ca_source_extraction'));
%     addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'utilities')));
%     addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'endoscope')));
%     addpath(fullfile(cnmfe_folder, 'GUI'));
%     addpath(fullfile(cnmfe_folder, 'GUI', 'gui_callbacks'));
%     addpath(fullfile(cnmfe_folder, 'GUI', 'modules'));
%     addpath(fullfile(cnmfe_folder, 'OASIS_matlab'));
%     addpath(fullfile(cnmfe_folder, 'scripts'));
%     cnmfe_loaded = true;
% end

%%
debug = 0;
if debug 
    visivility = 'on';
else
    visivility = 'off';
end

%% Select variables
animal_ID = 'MA_27';
sessions_to_take = [1, 2];
stimulus = 'SP';
% stimulus = 'HPS';
select_by_activity = false;
track_neurons = true;
take_raw = true; % Either raw or denoised
bkg_white = false; % Either background white or dark
%% Setting parameters
if strcmp(stimulus, 'SP')
    cell_idxs = [14, 6, 15, 5, 1]';
    move_n_frames = 15;
    vid_FR = 30;
    stimulus_at = [];
    trials_to_take = 5;
    what_trials = [1:trials_to_take];
else
    cell_idxs = [1,3:5,7:9,12:14]';
    move_n_frames = 3;
    vid_FR = 20;
    stimulus_at = 80;
    trials_to_take = 1;
    what_trials = 3;
end

if bkg_white
    line_color = 'b';
    bkg_color = 'w';
    str_disp = 'in white mode';
else
    line_color = 'w';
    bkg_color = 'k';
    str_disp = 'in dark mode';
end

fprintf('%% Making videos for %s %s %%\n',  stimulus, str_disp)
%% Load Traces
ROI_info_file = get_filename_of('ROI_info', animal_ID);
F = load_variable(ROI_info_file, 'ROI_fluorescence');

%% Load ROI info
disp('Loading ROIs')
ROI_path = os.path.join(GC.data_root_path, GC.segmentation_path);
ROI_info_file = os.path.join(ROI_path, [animal_ID, '_Neuron_Source2D.mat']);

%% Create Video Obj
rootpath = os.path.join(GC.temp_dir, animal_ID);
filename = os.path.join(rootpath,  ['FVO_video', '_', stimulus]);


%% load FVOs
disp('Loading Videos')
frames_idx =  get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',true);
METADATA = SQL_database.read_table_where('trials', {'+'}, animal_ID, 'animal_ID');
METADATA_sess = SQL_database.read_table_where('sessions', {'+'}, animal_ID, 'animal_ID');
days_str = (unique(METADATA_sess.day_from_surgery));
days_n = sort(str2double(days_str));
sessions_id = METADATA_sess.session_id(ismember(METADATA_sess.day_from_surgery, num2str(days_n(sessions_to_take))) & ismember(METADATA_sess.stimulus, stimulus));

% Get the frame ids to take
frames_idx_1 = find(ismember(METADATA.session_id, sessions_id(1)));
frames_idx_2 = find(ismember(METADATA.session_id, sessions_id(2)));
% select the initial frame and the last frame to read
% video 1
init_1= frames_idx(frames_idx_1,:);
init_1_idx = init_1(what_trials(1),1);
fin_1_idx = init_1(what_trials(end) ,2);
% video 2
init_2= frames_idx(frames_idx_2,:);
init_2_idx = init_2(what_trials(1),1); % +3 for HPS
fin_2_idx = init_2(what_trials(end),2);% +3 for HPS

% Load h5 file
h5_file = os.path.join(rootpath, [animal_ID, '_joint_corrected.h5']);
data_set = '/1';
h5_info = h5info(h5_file);
size_vid = h5_info.Datasets.Dataspace.Size(1:2);
n_frames  = fin_1_idx - init_1_idx +1;
% read file
FOV_1 = h5read(h5_file,data_set,[1,1,init_1_idx],[size_vid,n_frames]);

% Load h5 file
% read file
n_frames  = fin_2_idx - init_2_idx +1;
FOV_2 = h5read(h5_file,data_set,[1,1,init_2_idx],[size_vid,n_frames]);
disp('Videos loaded')
%% Get ROIs
neuron = load_variable(ROI_info_file, 'neuron_results');
disp('ROIs loaded')

%% Get traces
if take_raw
    TRACES = neuron.C_raw';
    do_downsample = 3;
else
    TRACES = neuron.C';
    do_downsample = 0;
end

%% Select neuron idxs
if select_by_activity % Will overwrite the cells alrady provided
    cell_idxs = [];
    act_neurons = mean(TRACES(init_1_idx+85:fin_1_idx,:),1)';
    [~, sorted_idx1 ] = sort(act_neurons, 'descend'); %#ok<*UDIM>
    % Session 2
    act_neurons = mean(TRACES( init_2_idx+85:fin_2_idx,:),1)';
    [~, sorted_idx2 ] = sort(act_neurons, 'descend');
    cell_idxs(:,:) = [sorted_idx1(1:5),sorted_idx2(1:5)];
    if ~track_neurons
        cell_idx1 = cell_idxs(:,1);
        cell_idx2 = cell_idxs(:,2);
    else
        cell_idx1 = cell_idxs(:,1);
        cell_idx2 = cell_idxs(:,1);
    end
else
    if ~track_neurons
        disp('Provide new neurons!!')
        keyboard
    else
        cell_idx1 = cell_idxs;
        cell_idx2 = cell_idxs;
    end
end

n_cells_to_show = size(cell_idxs,1);
% cells_to_take = cell_idxs;
Coor = neuron.Coor;
Aor1 = neuron.A(:,cell_idx1);
Aor2 = neuron.A(:,cell_idx2);
% Fluorecent traces
F_cells_1 = TRACES(init_1_idx:fin_1_idx, cell_idx1);
F_cells_2 = TRACES(init_2_idx:fin_2_idx, cell_idx2);

%% Create the figure
if ~isempty(stimulus_at)
    draw_square = 1;
    x = [0 round(size(FOV_1,1)/4) round(size(FOV_1,1)/4) 0];
    y = [0 0 round(size(FOV_1,1)/4) round(size(FOV_1,1)/4)];
else
    draw_square = 0;
end
close all
fig1 = figure('pos',[10 10 1500, 950], 'color', bkg_color, 'visible', visivility);
m = 2;n=2;
% Plot traces
% FOV1
ax3 = subplot(m,n,3);
ax3.Position = [0.13 0.125 0.4 0.37];
ylimits1 = 1500*ones(size(F_cells_1,1),1);
for i_cell = 1:n_cells_to_show
    data_this_cell = F_cells_1(:, i_cell) + ylimits1;
    if do_downsample > 0
        data_this_cell = smooth(data_this_cell', do_downsample)';
    end
    % Plot dFF
    hold on
    plot( data_this_cell * 100, 'Color',line_color, 'linewidth',.15);
    ylimits1 = ylimits1 + 3000;
    axis off
end
% FOV2
ax4 = subplot(m,n,4);
ax4.Position = [0.57 0.11 0.4 0.37];

ylimits2 = 1500*ones(size(F_cells_2,1),1);
for i_cell = 1:n_cells_to_show
    data_this_cell = F_cells_2(:, i_cell) + ylimits2;
    if do_downsample > 0
        data_this_cell = smooth(data_this_cell', do_downsample)';
    end
    % Plot dFF
    hold on
    plot( data_this_cell * 100, 'Color',line_color, 'linewidth',.15);
    ylimits2 = ylimits2 + 3000;
    axis off
end
%%
if debug, keyboard;end

%% Write video
tic
if ~debug
    writeObj = VideoWriter(filename, 'MPEG-4');
    writeObj.FrameRate = vid_FR;
    open(writeObj)
end
cmap_frame= 'gray';
iterations = round(linspace(1,n_frames, n_frames/move_n_frames));
if strcmp(visivility, 'off')
    f = waitbar(0,'Making video. Please wait...');
end
%%
for it = 1:length (iterations)-1
    %% FOV 1
    this_it = iterations(it);
    ax1 = subplot(m,n,1);
    plot_FOV_in_vid(mean(FOV_1(:,:,this_it:this_it+move_n_frames-1),3), Aor1, cmap_frame,cell_idx1, 1)
    if draw_square && this_it >= stimulus_at
        fill(ax1, x,y, 'r')
    end
    % FOV 2
    ax2 = subplot(m,n,2);
    plot_FOV_in_vid(mean(FOV_2(:,:,this_it:this_it+move_n_frames-1),3), Aor2, cmap_frame, cell_idx2, 2)
    if draw_square && this_it >= stimulus_at
        fill(ax2, x,y, 'r')
    end
    % TRACES 1
    li1 = plot(ax3, [this_it, this_it],[0, max(ylimits1)*100], 'r-');    
    %% TRACES 2
    li2 =plot(ax4, [this_it, this_it],[0, max(ylimits2)*100], 'r-');
    if ~debug
        % Get Frame
        Frame = getframe(fig1);
        % Write frame 
        writeVideo(writeObj, Frame)
    end
    % Delete lines
    delete(li1)
    delete(li2)
    if strcmp(visivility, 'off')
       waitbar(it/length (iterations), f,'Making video. Please wait...');
    end
end
time_elapsed = toc;
fprintf('It took %.2f mins \n', time_elapsed/60)
if strcmp(visivility, 'off')
    close(f)
end
close(writeObj)
close(fig1)
fprintf('Video is saved in %s \n', filename)

%% Helper Functions
% Function to draw ROIS
function plot_FOV_in_vid(data, Aor,cmap_frame, ~ , s)
imagesc(data), hold on
colormap(cmap_frame)
axis off
caxis([0 5000])
title(['Session ', num2str(s)])
d1 = size(data,1);
d2 = size(data,2);
thr = 0.995;
cmap = hot(3*size(Aor,2));
ln_wd = 2;
max_number = size(Aor,2);
fontname = 'Arial';
font_size = 16;
CC = cell(size(Aor,2),1);
CR = cell(size(Aor,2),2);
for i = 1:size(Aor,2)
    A_temp = full(reshape(Aor(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(size(Aor,2),:), 'linewidth', ln_wd);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end

% add text
cm = com(Aor(:,1:end),d1,d2);
lbl = strtrim(cellstr(num2str((1:size(Aor,2))')));
% lbl = cellstr(num2str(cell_idxs));
text((cm(1:max_number,2)),(cm(1:max_number,1)),lbl(1:max_number),'color','r','fontsize',font_size,'fontname',fontname,'fontweight','normal');
%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [srt] = orderROIs_vid(obj, srt)
% %% order neurons
% % srt: sorting order
% nA = sqrt(sum(obj.A.^2));
% nr = length(nA);
% if nargin<2
%     srt='srt';
% end
% K = size(obj.C, 1);
% 
% if ischar(srt)
%     if strcmpi(srt, 'decay_time')
%         % time constant
%         if K<0
%             disp('Are you kidding? You extracted 0 neurons!');
%             return;
%         else
%             taud = zeros(K, 1);
%             for m=1:K
%                 temp = ar2exp(obj.P.kernel_pars(m));
%                 taud(m) = temp(1);
%             end
%             [~, srt] = sort(taud);
%         end
%     elseif strcmp(srt, 'mean')
%         if obj.options.deconv_flag
%             temp = mean(obj.C,2)'.*sum(obj.A);
%         else
%             temp = mean(obj.C,2)'.*sum(obj.A)./obj.P.neuron.sn';
%         end
%         [~, srt] = sort(temp, 'descend');
%     elseif strcmp(srt, 'sparsity_spatial')
%         temp = sqrt(sum(obj.A.^2, 1))./sum(abs(obj.A), 1);
%         [~, srt] = sort(temp);
%     elseif strcmp(srt, 'sparsity_temporal')
%         temp = sqrt(sum(obj.C_raw.^2, 2))./sum(abs(obj.C_raw), 2);
%         [~, srt] = sort(temp, 'descend');
%     elseif strcmp(srt, 'circularity')
%         % order neurons based on its circularity
%         tmp_circularity = zeros(K,1);
%         for m=1:K
%             [w, r] = nnmf(obj.reshape(obj.A(:, m),2), 1);
%             ky = sum(w>max(w)*0.3);
%             kx = sum(r>max(r)*0.3);
%             tmp_circularity(m) = abs((kx-ky+0.5)/((kx+ky)^2));
%         end
%         [~, srt] = sort(tmp_circularity, 'ascend');
%     elseif strcmpi(srt, 'pnr')
%         pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);
%         [~, srt] = sort(pnrs, 'descend');
%     elseif strcmpi(srt, 'temporal_cluster')
%         obj.orderROIs('pnr');
%         dd = pdist(obj.C_raw, 'cosine');
%         tree = linkage(dd, 'complete');
%         srt = optimalleaforder(tree, dd);
%     elseif strcmpi(srt, 'spatial_cluster')
%         obj.orderROIs('pnr');
%         A_ = bsxfun(@times, obj.A, 1./sqrt(sum(obj.A.^2, 1)));
%         temp = 1-A_' * A_;
%         dd = temp(tril(true(size(temp)), -1));
%         dd = reshape(dd, 1, []);
%         tree = linkage(dd, 'complete');
%         srt = optimalleaforder(tree, dd);
%     elseif strcmpi(srt, 'snr')
%         snrs = var(obj.C_raw, 0, 2);%./var(obj.C-obj.C_raw, 0, 2);
%         [~, srt] = sort(snrs, 'descend');
%     end
% end
% obj.A = obj.A(:, srt);
% obj.C = obj.C(srt, :);
% 
% try
%     obj.C_raw = obj.C_raw(srt,:);
%     obj.S = obj.S(srt,:);
%     obj.P.kernel_pars = obj.P.kernel_pars(srt, :);
%     obj.P.neuron_sn = obj.P.neuron_sn(srt);
%     
%     obj.ids = obj.ids(srt);
%     obj.tags = obj.tags(srt);
% end
% end
