function ticklabels = log_ticklabels_fold_change(varargin)

% Parse inputs
p = inputParser();
addRequired(p, 'ticks', @(x) ~isempty(x))
addRequired(p, 'log_base', @(x) isnumeric(x) && isscalar(x) && (x > 0))
addOptional(p, 'return_type', 'fraction',  @(x) ismember(x, {'fraction', 'decimal', 'percent'}))
parse(p, varargin{:});
ticks = p.Results.ticks;
log_base = p.Results.log_base;
return_type = p.Results.return_type;

% Allocate output variable
n_ticks = length(ticks);
ticklabels = repmat({''}, n_ticks, 1);

% Convert positive ticks
positive_ticks = ticks > 0;
if any(positive_ticks)
    positive_ticks = find(positive_ticks);
    for i_tick = 1:length(positive_ticks)
        tick_value = ticks(positive_ticks(i_tick));
        switch return_type
            case 'fraction'
                tick_value_str = num2str(round(log_base ^ abs(tick_value)));
            case 'decimal'
                tick_value_str = num2str(log_base ^ abs(tick_value));
            case 'percent'
                tick_value_str = [num2str((log_base ^ abs(tick_value)) * 100), '%'];
        end
        ticklabels{positive_ticks(i_tick)} = tick_value_str;
    end
end

% Convert negative ticks
negative_ticks = ticks < 0;
if any(negative_ticks)
    negative_ticks = find(negative_ticks);
    for i_tick = 1:length(negative_ticks)
        tick_value = ticks(negative_ticks(i_tick));
        switch return_type
            case 'fraction'
                tick_value_str = ['1 / ', num2str(round(log_base ^ abs(tick_value)))];
            case 'decimal'
                tick_value_str = num2str(1 / (log_base ^ abs(tick_value)));
            case 'percent'
                tick_value_str = [num2str(1 / (log_base ^ abs(tick_value)) * 100), '%'];
        end
        ticklabels{negative_ticks(i_tick)} = tick_value_str;
    end
end

% Convert 0-tick
zero_tick = ticks == 0;
if any(zero_tick)
    switch return_type
        case 'fraction'
            tick_value_str = '1';
        case 'decimal'
            tick_value_str = '1';
        case 'percent'
            tick_value_str = '100%';
    end
    ticklabels{zero_tick} = tick_value_str;
end

