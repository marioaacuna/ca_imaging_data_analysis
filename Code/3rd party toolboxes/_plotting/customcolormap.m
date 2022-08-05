function cmap = customcolormap(stops, colors, n_colors, gradient_type)
    if ~exist('n_colors','var'), n_colors=256; end
    if ~exist('gradient_type','var'), gradient_type='lin'; end
    
    % Compute positions along the samples
    stops_position = round((n_colors-1) * stops) + 1;
    n_stops = length(stops_position);
    % Make the gradients among colors
    cmap = zeros(n_colors, 3);
    cmap(stops_position, :) = colors;

    n_intervals = n_stops - 1;
    for i_int = 1:n_intervals
        color_start = colors(i_int, :);
        color_end   = colors(i_int + 1, :);
        
        start_sample = stops_position(i_int);
        end_sample   = stops_position(i_int + 1);
        if start_sample < end_sample
            idx_range = start_sample:end_sample;
        else
            idx_range = start_sample:-1:end_sample;
        end
        n_shades = abs(end_sample - start_sample) + 1;
        
        for idx_rgb = 1:3
            switch gradient_type
                case 'lin'
                    color_gradient = linspace(color_start(idx_rgb), color_end(idx_rgb), n_shades);
                case 'log'
                    color_gradient = logspace(log10(color_start(idx_rgb)), log10(color_end(idx_rgb)), n_shades);
            end
            cmap(idx_range, idx_rgb) = color_gradient;
        end
    end

    cmap = cmap ./ 255;
