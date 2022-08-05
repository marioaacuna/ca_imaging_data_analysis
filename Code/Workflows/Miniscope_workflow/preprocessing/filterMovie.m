function movie = filterMovie(p, movie, fType, remove_nans)
%FILTERMOVIE filters movie according to the paramneters set in p

if ~exist('remove_nans', 'var'), remove_nans = true; end

% remove nans
if remove_nans
    movie(isnan(movie)) = 0;
end

%% create cutoffFilter and normalize to be between 0 and 1
padSize = size(movie,1);
filterSize = [size(movie,1) + 2*padSize , size(movie,2) + 2*padSize];

switch fType 
    case 'lowpass'
        cutoffFilter = mat2gray(fspecial('gaussian', filterSize, p.filtering.lowpassFreq));

    case 'highpass'
        cutoffFilter = 1 - mat2gray(fspecial('gaussian', filterSize, p.highpassFreq));

    case 'bandpass'
        lowpassFilter = mat2gray(fspecial('gaussian', filterSize, p.turboreg.bandpassFreqs(2)));
        highpassFilter = 1 - mat2gray(fspecial('gaussian', filterSize, p.turboreg.bandpassFreqs(1)));
        cutoffFilter = highpassFilter .* lowpassFilter;
end

%% filter all frames
% overwrite in loop to reduce memory
parfor i = 1:size(movie,3)
    movie(:, :, i) = filterImage(movie(:, :, i), cutoffFilter, padSize, fType);
end
