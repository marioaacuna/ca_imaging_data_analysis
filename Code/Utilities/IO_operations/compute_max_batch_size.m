function batch_size = compute_max_batch_size(start_value, max_percent_memory, tested_on)

batch_size = start_value;
% Adapt batch size to current memory capacity to avoid out-of-memory errors
if ispc()
    % Get memory info
    [~, systemview] = memory;
    current_RAM = systemview.PhysicalMemory.Available / 1024 / 1024 / 1024;
    available_RAM = (max_percent_memory/100) * current_RAM;  % Do not overload memory so Windows will not freeze
    batch_size = floor(tested_on * available_RAM / 32);
end
