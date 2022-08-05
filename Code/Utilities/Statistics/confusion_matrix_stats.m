function stats = confusion_matrix_stats(confusion_matrix)

n_classes = size(confusion_matrix, 1);

accuracy = NaN(n_classes, 1);
sensitivity = NaN(n_classes, 1);
specificity = NaN(n_classes, 1);
precision = NaN(n_classes,1);
matthews_correlationcoefficient =  NaN(n_classes,1);
F1 = NaN(n_classes, 1);

for i_class = 1:n_classes
    is_class = i_class;
    is_not_class = setdiff(1:n_classes, is_class);
    
    % Calculate proportions
    n = nansum(confusion_matrix(:));
    TP = confusion_matrix(is_class, is_class);
    TN = confusion_matrix(is_not_class, is_not_class); 
    TN = nansum(TN(:));
    FN = nansum(confusion_matrix(is_class, is_not_class));
    FP = nansum(confusion_matrix(is_not_class, is_class));
    
    accuracy(i_class) = (TP + TN) / n;
    sensitivity(i_class) = TP / (TP + FN);
    specificity(i_class) = TN / (TN + FP);
    F1(i_class) = (2 * TP) / (2 * TP + FP + FN);
    precision(i_class) = TP / (TP + FP);
    matthews_correlationcoefficient(i_class) = (TP * TN - FP * FN) / (sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)));
    
end

% Gather results
stats = struct();
stats.accuracy = accuracy;
stats.sensitivity = sensitivity;
stats.specificity = specificity;
stats.F1 = F1;
stats.precision = precision;
stats.MCC = matthews_correlationcoefficient;
stats.avg_accuracy = nanmean(accuracy);
stats.avg_sensitivity = nanmean(sensitivity);
stats.avg_specificity = nanmean(specificity);
stats.avg_F1 = nanmean(F1);
stats.avg_precision = nanmean(precision);
stats.avg_MCC = nanmean(matthews_correlationcoefficient);

