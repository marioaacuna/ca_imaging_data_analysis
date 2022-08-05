function varargout = make_CDF(values, bins, return_PDF)
    PDF = hist(values, bins);
    PDF = PDF ./ length(values);
    % If requested, return PDF
    if return_PDF
        varargout{3} = [PDF(:), bins(:)];
    end
    
    % Compute the cumulative distribution function (CDF)
    CDF = cumsum(PDF);
    % Keep only unique bins
    unique_idx = [1, find(diff(CDF)>0)+1];
    CDF = CDF(unique_idx);
    bins = bins(unique_idx);
    % Assign outputs
    varargout{1} = CDF;
    varargout{2} = bins;

