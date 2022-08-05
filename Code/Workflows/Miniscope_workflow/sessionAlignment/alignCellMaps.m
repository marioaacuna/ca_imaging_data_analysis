function [ mapsAligned, tforms, stack ] = alignCellMaps(p, maps)
% generates tforms (transformation matrix) for each session s.t. all
% sessions from one mouse can be aligned. this is a bit tricky because 
% turboreg has to be run iteratively, generating multiple sets of 
% registration coords. These coords can't be accurately combined to one set,
% therefore we additionaly use a matlab registration algorithm which 
% takes the (almost correct) registration coords generated by turboreg and 
% does the fine tuning s.t. we get one tform that can be used to align the 
% session movies.


%% Align maps with turboreg
mapsTR = cat(3, maps{:});

pAlign = p;
pAlign.turboreg.refFrame = 1; 
pAlign.turboreg.options.RegisType = 2;
pAlign.turboreg.options.Levels = 6;

nIter = 5;
regCoords = cell(1, nIter);
for i_iter = 1:nIter
    [mapsTR, regCoords{1, i_iter}] = turboregCellMaps(pAlign, mapsTR);
end


%% Get final registration coordinates from iterative alignment
getTransMat = @(coords) [1 0 coords(1); 0 1 coords(2); 0 0 1];
getRotMat = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

n_sessions = length(maps);
tformsTR = cell(1, n_sessions);
for i_sess = 1:n_sessions
    % reset pMat for every session
    pMat = eye(3);
    for i_iter = 1:nIter
        % get values from this iteration
        trans = flipud(regCoords{i_iter}{i_sess}.Translation);
        ori = flipud(regCoords{i_iter}{i_sess}.Origin);
        rot = regCoords{i_iter}{i_sess}.Rotation;      
        % calculate projection matrix 
        pMat = pMat * getTransMat(ori) * getRotMat(rot) * getTransMat(-ori) * getTransMat(trans);
    end
    tformsTR{i_sess} = affine2d(pMat'); 
end


%% Align with imreg, using truboreg tforms as inital transformation
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
optimizer.InitialRadius = optimizer.InitialRadius/3.5;

mapsAligned = cell(1, n_sessions);
tforms = cell(1, n_sessions);
fixed = maps{1}; % align to day 1
for i_sess = 1:n_sessions
    moving = maps{i_sess};
    tforms{i_sess} = imregtform(moving, fixed, 'rigid', optimizer, metric,'InitialTransformation', tformsTR{i_sess});
    mapsAligned{i_sess} = imwarp(maps{i_sess}, tforms{i_sess}, 'OutputView',imref2d(size(maps{i_sess})));
end


%% Inspect all aligned maps 
mapSizes = cell2mat(cellfun(@(x) size(x), mapsAligned, 'UniformOutput',false));
mapSizes = reshape(mapSizes, [2, n_sessions])';
maxX = max(mapSizes(:, 1));
maxY = max(mapSizes(:, 2));

stack = zeros(maxX, maxY, n_sessions, 'single');
for i_sess = 1:n_sessions
    movieSize = size(mapsAligned{i_sess});
    stack(1:movieSize(1), 1:movieSize(2), i_sess) = mapsAligned{i_sess};
end
imtool3D(stack);


end

