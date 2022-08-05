function [thisMovie] = dfofMovie(thisMovie)
% converts movie to units of DFOF
    
    % adjust for problems with movies that have negative pixel values before dfof
    minMovie = min(thisMovie(:));
    if minMovie < 0
        thisMovie = thisMovie + 1.1*abs(minMovie);
    end

    % get the movie F0, do by row to reduce potential memory errors	    
    inputMovieF0 = zeros([size(thisMovie,1) size(thisMovie,2)]);
    parfor rowNo = 1:size(thisMovie,1)
        inputMovieF0(rowNo,:) = nanmean(squeeze(thisMovie(rowNo,:,:)),2);
    end
    
    % bsxfun for fast matrix divide
    thisMovie = bsxfun(@ldivide,inputMovieF0,thisMovie);
    thisMovie = thisMovie-1;

end


