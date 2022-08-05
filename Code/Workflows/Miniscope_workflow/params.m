function p = params()
% PARAMS builds structure containing all parameters 
% all paramters have a default value, some should be adjusted to the data 


%% filtering
    p.filtering.lowpassFreq = 4;  % see preprocessing/testFilters for exploring different values

%% turboreg
    p.turboreg.options.RegisType=1;  % 1 is only translation, see turboreg website for other types
    p.turboreg.options.SmoothX=10;
    p.turboreg.options.SmoothY=10;
    p.turboreg.options.minGain=0.0;
    p.turboreg.options.Levels=4; 
    p.turboreg.options.Lastlevels=1;
    p.turboreg.options.Epsilon=1.192092896E-07;
    p.turboreg.options.zapMean=0;
    p.turboreg.options.Interp='bilinear';
    
    p.turboreg.refFrame = 100;
    p.turboreg.bandpassFreqs = [20 80];  % see preprocessing/testFilters for exploring different values
    
%% downsample time
    p.downsampleTime.secondaryDownsampleType = 'bilinear';
    
%% PCA ICA
    p.PCAICA.mu = 0.1;
    p.PCAICA.term_tol = 1e-5;
    p.PCAICA.max_iter = 750;
    
end

