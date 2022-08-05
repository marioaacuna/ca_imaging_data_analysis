function [data_class0, data_class1] = standardize_data(data_class0, data_class1)
% Get values for normalization
data = cat(2, data_class1, data_class0);
mean_activity = mean(data, 2);
std_activity = std(data, 0, 2);
% Apply normalization
data_class1 = (data_class1 - mean_activity) ./ std_activity;
data_class0 = (data_class0 - mean_activity) ./ std_activity;
% Replace possible NaNs with 0s
data_class1(~isfinite(data_class1)) = 0;
data_class0(~isfinite(data_class0)) = 0;
