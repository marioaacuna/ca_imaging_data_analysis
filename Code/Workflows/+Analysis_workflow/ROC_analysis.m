function varargout = ROC_analysis(data, varargin)
% Function "roc_curve" calculates the Receiver Operating Characteristic  
% curve, which represents the 1-specificity and sensitivity, of two classes
% of data, called class_1 and class_2.                                    
%                                                                         
%   Input parameters                                                      
%       · data:     Two-column matrix which stores the data of both classes
%                   with the following structure:                          
%                   data = [class_1, class_2]                           
%                   where class_1 and class_2 are column vectors.         
%       · dispp:    (Optional) If dispp is 1, the ROC Curve will be displayed
%                   inside the active figure. If dispp is 0, no figure
%                   will be displayed.                                    
%       · dispt:    (Optional) If dispt is 1, the optimum threshold para
% 					meters obtained will be displayed on the MATLAB log.  
%                   Otherwise, if dispt is 0, no parameters will be displayed
% 					there.                                           
%                                                                         
%   Output variables                                                      
%       · ROC_data: Structure that contains all the curve parameters.     
%           - results:  Structure that contains the quantitative parameters 
%                   of the obtained curve, which are:                 
%               + Threshold: Optimum threshold calculated in order to maximize
% 					the sensitivity and specificity values, which is colocated
% 					in the nearest point to (0,1).                                        
%               + AROC:     Area Under ROC Curve.                         
%               + Accuracy: Maximum accuracy obtained.                    
%               + Sensitivity:    Optimum threshold sensitivity.                
%               + Specificity:    Optimum threshold specificity.                
%               + PPV:      Positive predicted value.                     
%               + NPV:      Negative predicted value.                     
%           - curve:    Matrix which contains the specificity and specifi-
%                       city of each threshold point in columns.          
% ----------------------------------------------------------------------- 
%   Example use:                                                         
%       class_1 = 0.5*randn(100,1);                                       
%       class_2 = 0.5+0.5*randn(100,1);                                   
%       data = [class_1, class_2];                                        
%       ROC_Curve(data);                                                  

    % Parse inputs
    p = inputParser();
    addParameter(p, 'class_order_matters', true)
    addParameter(p, 'show_figure', false)
    addParameter(p, 'print_output', false)
    parse(p, varargin{:})
    % Unpack paramaters
    class_order_matters = p.Results.class_order_matters;
    show_figure = p.Results.show_figure;
    print_output = p.Results.print_output;
    
    % Get data dimension
    L = size(data, 1);
 
    % Allocate output variable
    results = struct();
    results.threshold = NaN;
    results.sensitivity_thr = NaN;
    results.specificity_thr = NaN;
    results.TP = NaN;
    results.TN = NaN;
    results.FP = NaN;
    results.FN = NaN;
    results.TPR = NaN;
    results.SPC = NaN;
    results.PPV = NaN;
    results.NPV = NaN;
    results.FPR = NaN;
    results.FNR = NaN;
    results.FDR = NaN;
    results.accuracy = NaN;
    results.F1_score = NaN;
    results.MCC = NaN;
    results.informedness = NaN;
    results.markedness = NaN;
    results.AUROC = 0.5;
    results.G = NaN;
    
    % Calculate the threshold values between data points
    s_data = unique(sort(data(:)));
    d_data = diff(s_data);
    if isempty(d_data)
        if nargout >= 1
            varargout{1} = results;
        end
        return
    end
    % Compute range of thresholds
    d_data(length(d_data)+1,1) = d_data(length(d_data));
    thres = NaN(size(s_data,1), 1);
    thres(1) = s_data(1) - d_data(1);
    thres(2:length(s_data)+1) = s_data + d_data ./ 2;

    % Sort classes
    if ~class_order_matters
        if mean(data(:,1)) > mean(data(:,2))
            data = [data(:,2), data(:,1)];
        end
    end
    
    % Calculating the sensibility and specificity at each threshold
    curve = zeros(size(thres,1),2);
    distance = zeros(size(thres,1),1);
    for id_t = 1:1:length(thres)
        TP = length(find(data(:,2) >= thres(id_t)));    % True positives
        FP = length(find(data(:,1) >= thres(id_t)));    % False positives
        FN = L - TP;                                    % False negatives
        TN = L - FP;                                    % True negatives
        
        curve(id_t,1) = TP/(TP + FN);   % Sensitivity
        curve(id_t,2) = TN/(TN + FP);	% Specificity
        
        % Distance between each point and the optimum point (0,1)
        distance(id_t) = sqrt((1-curve(id_t,1)).^2 + (curve(id_t,2)-1).^2);
    end
    
    % Optimum threshold and parameters
    [~, opt] = min(distance);
    % Recalculate values of confusion matrix
    TP = length(find(data(:,2) >= thres(opt)));
    FP = length(find(data(:,1) >= thres(opt)));
    FN = L - TP;
    TN = L - FP;
    
    % Store outputs
    % Optimum threshold position
    results.threshold = thres(opt);
    % Optimum threshold's sensitivity and specificity
    results.sensitivity_thr = curve(opt,1);
    results.specificity_thr = curve(opt,2);
    % Store confusion matrix
    results.TP = TP;
    results.TN = TN;
    results.FP = FP;
    results.FN = FN;
    % Derive performance from confusion matrix (https://en.wikipedia.org/wiki/Sensitivity_and_specificity#Definitions)
    results.TPR = TP/(TP+FN);  % sensitivity
    results.SPC = TN/(TN+FP);  % specificity
    results.PPV = TP/(TP+FP);  % positive predictive value or precision
    results.NPV = TN/(TN+FN);  % negative predictive value
    results.FPR = FP/(FP+TN);  % false positive rate or fall-out
    results.FNR = FN/(TP+FN);  % false negative rate
    results.FDR = FP/(TP+FP);  % false discovery rate
    results.accuracy = (TP+TN)/(TP+FP+FN+TN);
    results.F1_score = 2*TP/(2*TP+FP+FN);  % harmonic mean of precision and sensitivity
    results.MCC = (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));  % Matthews correlation coefficient between the observed and predicted binary classifications (also known as phi coefficient)
    results.informedness = results.TPR + results.SPC - 1;  % Youden's J statistic
    results.markedness = results.PPV + results.NPV - 1;
    % Area-under-the-ROC curve
    results.AUROC = abs(trapz(1-curve(:,2), curve(:,1)));
    % Gini's coefficient
    results.G = 2 * results.AUROC - 1;
        
    %% Plot
    if show_figure
        figure('color','w')
        clf, hold on
        plot([0,1], [0,1], 'color',[.7,.7,.7], 'LineWidth',2)
        plot(1-curve(:,2), curve(:,1), 'b', 'LineWidth',2)
        plot(1-curve(opt,2), curve(opt,1), 'or', 'MarkerSize',10, 'Linewidth',3)
        set(gca, 'Layer','top', 'TickDir','out', 'XTick',0:.1:1, 'YTick',0:.1:1)
        axis equal square; grid on
        xlabel('False Positive Rate'), ylabel('True Positive Rate')
        title(['AUROC = ' num2str(results.AUROC)])
    end
    
    %% Print outputs
    if print_output
        disp(results)
    end
    
    %% Return outputs
    if nargout >= 1
        varargout{1} = results;
    end
    if nargout >= 2
        varargout{2} = curve;
    end
end
