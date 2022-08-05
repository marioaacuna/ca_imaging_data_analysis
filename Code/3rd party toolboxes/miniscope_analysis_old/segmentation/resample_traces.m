function traces_out = resample_traces(traces_in, upsampling_factor, downsampling_factor)

traces_in = double(traces_in);

n_cells = size(traces_in, 1);

y = resample(traces_in(1, :), upsampling_factor, downsampling_factor);
new_n_frames = length(y);
traces_out = NaN(n_cells, new_n_frames);

for i_cell = 1:n_cells
    traces_out(i_cell, :) = resample(traces_in(i_cell, :), upsampling_factor, downsampling_factor);
end
     