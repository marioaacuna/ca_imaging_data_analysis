function [ movie,zeroMask ] = registerMovie(p, movie, verbose)
%REGISTERMOVIE registeres movie using turboreg. first registration 
%coordinates are obtained on the a bandpass filtered and inverted version
%of the movie. The original movie is then registered with these coordinates  

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

if verbose
    global LOGGER
end

%% preprocess movie to improve turboreg performance
if verbose
    LOGGER.info([char(datetime('now')) ' registering movie: starting filtering'])
end

% bandpass filter
movieFiltered = filterMovie(p,movie,'bandpass');

% invert, normalize and mean subtract frames
parfor i = 1:size(movieFiltered,3)
    movieFiltered(:,:,i) = mat2gray(imcomplement(movieFiltered(:,:,i)));
    frame = squeeze(movieFiltered(:,:,i));
    frameMean = mean(frame(:));
    movieFiltered(:,:,i) = frame - frameMean;
end
        
%% Get registration coordinates with turboreg
if verbose
    LOGGER.info([char(datetime('now')) ' registering movie: starting turboreg'])
end

% set params 
options = p.turboreg.options;

% get reference frame and set mask (just a dummy variable that is not
% really used)
refFrame = p.turboreg.refFrame;
if refFrame > size(movieFiltered, 3)
    refFrame = size(movieFiltered, 3);
end
refPic = single(movieFiltered(:, :, refFrame));
%refPic = single(movieFiltered);
mask = zeros(size(refPic),'single'); 

% run turboreg
parfor i = 1:size(movieFiltered,3)
    pic = single(movieFiltered(:,:,i));
    [~, regCoords{i}] = turboreg(refPic,pic,mask,mask,options);
end
    
    
%% Perform registration with transfturboreg
if verbose
    LOGGER.info([char(datetime('now')) ' registering movie: starting registration'])  
end
    
mask = ones(size(movie(:,:,1)),'single');
parfor i = 1:size(movie,3)
    pic = single(movie(:,:,i));
    movie(:,:,i) = transfturboreg(pic,mask,regCoords{i});
end    
    
%% get mask of zero pixels caused by translation 
% used for cropping at the end of preprocessing

zeroMask = squeeze(max(movie == 0,[],3));


%#ok<*TLEV>
