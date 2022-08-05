function values = compute_ICDF(N, CDF_y, CDF_x, interpolate)
    
    % Make sure to have column vectors
    CDF_y = CDF_y(:);
    CDF_x = CDF_x(:);

    % Draw random valeus from the continuous uniform distribution
    r = rand(N,1);

    if length(CDF_y)==1 && length(CDF_x)==1
        % Simply repeat the only known value
        values = repmat(CDF_y, N, 1);
        
    else
        if ~interpolate  % Simply get the value that is closest in the CDF
            [~, idx] = min(abs(CDF_y(:)' - r), [] ,2);
            values = CDF_x(idx);
            
            % Make sure output is column vector
            values = values(:);
        else
            try
                % Keep only unique bins
                unique_idx = [1; find(diff(CDF_y)>0)+1];
                CDF_y = CDF_y(unique_idx);
                CDF_x = CDF_x(unique_idx);
                % Find values by lienar interpolatation
                values = interp1(CDF_y, CDF_x, r);
                % Replace NaNs with the lowest value
                values(isnan(values)) = CDF_x(1);
            catch ME
                if strcmpi(ME.identifier, 'MATLAB:griddedInterpolant:DegenerateGridErrId')
                    values = zeros(N,1);
                end
            end 
        end
    end
