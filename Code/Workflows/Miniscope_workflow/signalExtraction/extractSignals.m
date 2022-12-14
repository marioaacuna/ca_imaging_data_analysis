function [ output_args ] = extractSignals( p, moviePath,saveDir )
%EXTRACTSIGNALS

    % load and crop movie
    movie = loadMovie(moviePath);
    %movie = load (moviePath, 'I_sq');
    zeroMask = squeeze(max(movie == 0,[],3));
    [N, S, W, E] = getCropCoords(zeroMask);
    movie = movie(N:S,W:E,:);   
    
    % run PCAICA
    [filters, traces] = runPCAICA(p, movie);
   
    % put filters back into right frame and save
    fFrame = zeros([size(zeroMask) size(filters,3)]);
    fFrame(N:S,W:E,:) = filters;
    filters = fFrame;
    cellMap = squeeze(max(filters,[],3));

    % save data
    save( [saveDir '\resultsPCAICA'],'p','traces','filters');
    save( [saveDir '\cellMap'],'cellMap');

end

