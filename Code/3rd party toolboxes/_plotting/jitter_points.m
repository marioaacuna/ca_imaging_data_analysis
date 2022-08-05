function jitter_vector = jitter_points(n_points, max_jitter)

% Run demo
if nargin == 0, run_demo(), return, end

%% FUNCTION
% Solve the trivial case
if n_points == 1
    jitter_vector = 0;
    return
end

% Compute the total distance over which we are allowed to space the points
total_distance = max_jitter * 2;

% The number of segments in which to divide the total distance equals the number
% of points to place at the center of each segment
n_segments = n_points;
    
% Compute the length of each segment
segment_length = total_distance / n_segments;

% Because each point is located in the middle of its segment, the extrema are
% not located anymore at the maximum distance, but they are half segment_length
% closer to 0
left_extremum  = -max_jitter + segment_length / 2;
right_extremum = +max_jitter - segment_length / 2;

% Compute linearly spaced positions
jitter_vector = linspace(left_extremum, right_extremum, n_points);


%% DEMO
function run_demo()
    % Open figure and axes
    figure('color','w')
    hold on
    
    % Set the maximum jitter to add to each point in both directions
    max_jitter = 1;
    
    % Loop and show how points are jittered accordingly to their number
    for n_points = 1:10
        j = jitter_points(n_points, max_jitter);
        plot(zeros(n_points, 1) + j, ones(n_points, 1) .* n_points, 'ok', 'markerfacecolor','k', 'markersize',14)
    end
    
    % Fix axes appearance
    xlabel('Point position')
    ylabel('Number of points')
    title('Demostration of jitter_points()', 'Interpreter','none')
    set(gca, 'XLim',[-max_jitter, max_jitter], 'YGrid','on', 'YLim',[0, 11], 'YTick',1:10, 'TickDir','out')

