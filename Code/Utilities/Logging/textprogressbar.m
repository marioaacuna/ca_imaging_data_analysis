function textprogressbar(c, n_iter)
% This function creates a text progress bar. It should be called with a 
% STRING argument to initialize and terminate. Otherwise the number correspoding 
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           % Carriage return pesistent variable
persistent nIter;           % Number of iterations to reach 100%
persistent lastIteration;   % Progress of last iteration

% Vizualization parameters
strPercentageLength = 7;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main 

if ~exist('c','var') || isempty(c)
    c = lastIteration + 1;
end

% Progress bar - initialization
if ischar(c) && isnumeric(n_iter) || isempty(strCR)
    if ~ischar(c)
        fprintf('progress: ');
    else
        fprintf('%s ',c);
    end
    strCR = -1;
    nIter = n_iter;
    lastIteration = 0;
    
% Progress bar  - termination
elseif ~isempty(strCR) && ischar(c)
    strCR = [];
    nIter = 1;
    fprintf([c '\n']);
    
% Progress bar - normal progress
elseif isnumeric(c)
    lastIteration = c;
    is_last_iteration = c >= nIter;
    c = floor(100 / nIter * c);
    if c < 0, c = 0; end
    if c > 100, c = 100; end
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
    % Automatic termination
    if is_last_iteration
        strCR = [];
        nIter = 1;
        lastIteration = 0;
        fprintf('\n');
    end
    
else
    % Any other unexpected input
    error('progress:unknown_error', 'Unsupported argument type');
end
