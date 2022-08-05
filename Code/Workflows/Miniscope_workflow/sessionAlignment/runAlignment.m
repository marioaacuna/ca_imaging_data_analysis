function runAlignment(animal_ID, experiments)
%runAlignment aligns cell maps of all sessions for each subject and saves
%registration tform

global GC LOGGER


% Load parameters
p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
p = load_variable(p_filename, 'p');

root_dir_on_server = os.path.join(GC.registered_images, p.subjects{1});
filename = os.path.join(root_dir_on_server, [p.subjects{1}, '_alignedCellMaps.mat']);

if exist(filename, 'file')
    LOGGER.info('Alignment done, skiping')
    return
end

LOGGER.info('Aligning cell maps of all sessions')

%% Load maps
n_sessions = size(experiments, 1);
maps = cell(1, n_sessions);
for i_sess = 1:n_sessions
    this_experiment = experiments(i_sess, :);
    
    filename = get_filename_of('miniscope_cell_map', p.subjects{1}, this_experiment{1}, this_experiment{2});
    maps{i_sess} = load_variable(filename, 'cellMap');
end


%% Run alignment
[mapsAligned, tforms, stack] = alignCellMaps(p, maps);


%% Save tforms and aligned maps
save(filename, 'mapsAligned', 'stack', 'tforms');


%#ok<*ASGLU>
