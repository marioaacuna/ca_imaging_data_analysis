function [p, h, stats] = wilcoxon_signrank(x,y,varargin)

try
    [p, h, stats] = signrank(x, y, varargin{:});

catch ME
    if strcmpi(ME.identifier, 'stats:signrank:NotEnoughData')
        p = NaN;
        h = NaN;
        stats = struct();
        
    else
        rethrow(ME)
    end
end
        

