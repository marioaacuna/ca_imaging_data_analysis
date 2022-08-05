function CCimage = CrossCorrImage(data, w, n_frames_per_CC_image)

if ~exist('w','var') || isempty(w)
    w = 1; % Window size
end

% Get size of video
video_height = size(data,1);
video_width = size(data,2);
tot_n_frames = size(data,3);

if ~exist('n_frames_per_CC_image','var') || isempty(n_frames_per_CC_image)
    n_frames_per_CC_image = tot_n_frames;
end

% Get number of chunks
n_chunks = ceil(tot_n_frames / n_frames_per_CC_image);

% Allocate output matrix
CCimage = zeros(video_height, video_width, n_chunks);

% Loop through chunks
for ichunk = 1:n_chunks
    % Get edges of chunk (and account for last chunk)
    start_sample = (ichunk-1) * n_frames_per_CC_image + 1;
    end_sample = min([start_sample + n_frames_per_CC_image, tot_n_frames]);
    % Slice data
    this_data = im2double(data(:, :, start_sample:end_sample));
    n_frames_in_chunk = size(this_data, 3);
    
    % Allocate matrix for computation
    ccimage = zeros(video_height, video_width);
    
    % Compute crosscorrelation
    for y = 1+w : video_height-w
        for x = 1+w : video_width-w
            % Center pixel
            a = this_data(y,x,:);
            center_pixel = reshape(a - mean(a, 3), [1, 1, n_frames_in_chunk]);  % Extract center pixel's time course and subtract its mean
            ac_center_pixel = sum(center_pixel.^2, 3);  % Autocorrelation, for normalization later

            % Neighborhood
            a = this_data(y-w:y+w, x-w:x+w, :);  % Extract the neighborhood
            neighboring_pixels = bsxfun(@minus, a, mean(a, 3));  % Subtract its mean
            ac_neighboring_pixels = sum(neighboring_pixels.^2, 3);  % Autocorrelation, for normalization later

            % Crosscorrelation
            product_pixels = bsxfun(@times, center_pixel, neighboring_pixels);
            product_ac = bsxfun(@times, ac_center_pixel, ac_neighboring_pixels);

            ccs = sum(product_pixels, 3) ./ sqrt(product_ac); % Crosscorrelation with normalization
            ccs((numel(ccs)+1) ./ 2) = [];  % Delete the middle point
            ccimage(y,x) = mean(ccs(:));  % Get the mean cross-correlation of the local neighborhood
        end
    end
    
    % Store this matrix
    CCimage(:, :, ichunk) = ccimage;
end
