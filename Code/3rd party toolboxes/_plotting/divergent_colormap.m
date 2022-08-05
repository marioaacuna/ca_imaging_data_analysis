function varargout = divergent_colormap (value_min, value_zero, value_max, monochrome, invert, cmapFunc, return_n_zero)

if ~exist('monochrome','var')
    monochrome = false;
end
if ~exist('cmapFunc','var')
    cmapFunc = [];
end
if ~exist('invert','var')
    invert = true;
end
if ~exist('return_n_zero','var')
    return_n_zero = false;
end

if isempty(cmapFunc)
    if monochrome
        cmapFunc = @(x) colormap(gray(x));
    else
        cmapFunc = @(x) cbrewer('div','RdBu',x);
    end
end    
    
MaxAllowedColorShades = 256;
color_step = (value_max-value_min) / MaxAllowedColorShades;

PositiveColors = ceil(abs(value_max-value_zero) / color_step);
NegativeColors = MaxAllowedColorShades - PositiveColors - 1;
halfSpectrum = max([PositiveColors, NegativeColors]);
Colors2rem = abs(PositiveColors-NegativeColors);
    
CMAP = cmapFunc(halfSpectrum*2+1);
if invert, CMAP = flipud(CMAP); end

% Remove out-of-bounds colors
if PositiveColors > NegativeColors
    CMAP(1:Colors2rem,:) = [];
else
    CMAP(end-Colors2rem+1:end,:) = [];
end

varargout{1} = CMAP;
if return_n_zero
    varargout{2} = NegativeColors;
end
