function test_FOV_correction(animal_ID, exp_type)
% test FVO corrected file
% animal_ID = 'MA_35epi_a';
filename = ['D:\_MATLAB_CaImaging\', animal_ID, '\', animal_ID, '_joint_corrected.h5'];
h5_info = h5info(filename);
n_frames = h5_info.Datasets.Dataspace.Size(3);
data_size = h5_info.Datasets.Dataspace.Size;
all_frames = [];
take_precent = 35; % in percent %
n_frames_to_take = round((n_frames-200) / take_precent);
frames_to_take = round(linspace(1,n_frames, n_frames_to_take));
figure('pos', [10 10 1200 900])
ex_frame_idx = frames_to_take(100);
ex_frame = h5read(filename, '/1' ,[1,1, ex_frame_idx], [data_size(1:2), 1]);
max_ax = max(max(ex_frame));
min_ax = min(min(ex_frame));
if strcmp(exp_type, 'epi')
    caxis_val = [min_ax,max_ax];
else
   caxis_val = [min_ax,max_ax];
end
for i_f = 1:n_frames_to_take-4
    this_frame = frames_to_take(i_f);
    this_whole_frames =  h5read(filename, '/1' ,[1,1, this_frame], [data_size(1:2), 100]);
    this_whole_frame = mean(this_whole_frames,3);
    all_frames(:,:,i_f) = this_whole_frame;
    imagesc(this_whole_frame) , caxis([caxis_val])
    title(num2str(this_frame))
    drawnow
end
end
%#ok<*AGROW,*SAGROW>
% implay(all_frames)