# Compute change-points, if requested
if PARAMETERS['classification_change_point_detection']:
    temp_dir = tempfile.gettempdir()
    R_input_filename = os.path.join(temp_dir, 'cpm.csv')
    R_log = os.path.join(temp_dir, 'R.log')
    current_folder = os.path.dirname(os.path.realpath(__file__))
    R_script = os.path.join(current_folder, 'decoding', 'detect_changepoint.R')

RESULTS['class_comparisons']['change_points'] = dict()
for class_idx, cl in enumerate(class_names):
    if cl == 'SP':
        res = []

    else:
        try:
            # Get indices of "true positive" assignments
            TP = np.where((SAMPLE_IDX['class'] == class_idx) & (
                        SAMPLE_IDX['predicted_label'].values.astype(
                            int) == class_idx) & (SAMPLE_IDX[
                                                      'significant_classified'] == 1))[
                0]
            sample_idx_in_class = SAMPLE_IDX.loc[
                TP, 'sample_idx_in_class'].values.astype(int)
            # Get data in those bins
            activity = np.mean(PARAMETERS['data'][cl], axis=1).reshape(
                    n_trials_per_class[class_idx], -1).ravel()
            y = np.zeros_like(activity) * np.nan
            y[sample_idx_in_class] = activity[sample_idx_in_class]
            y = y.reshape(n_trials_per_class[class_idx], -1)
            y = np.nanmean(y, axis=0)
            y[np.isnan(y)] = 0

            # ACF and PACF
            from statsmodels.tsa.stattools import acf, pacf

            lag_acf = acf(y, nlags=y.shape[0])
            lag_pacf = pacf(y, nlags=y.shape[0] - 1, method='ols')

            # Number of AR (Auto-Regressive) terms (p): AR terms are just lags of
            # dependent variable. For instance if p is 5, the predictors for x(t)
            # will be x(t-1)….x(t-5).
            # p – The lag value where the PACF plot crosses the upper confidence
            # interval for the first time.
            uci = 1.96 / np.sqrt(len(y))
            p = np.where(lag_pacf <= uci)[0]
            if p.size > 0:
                p = p[0] + 1
            else:
                p = 1

            # Number of MA (Moving Average) terms (q): MA terms are lagged forecast
            # errors in prediction equation. For instance if q is 5, the predictors
            # for x(t) will be e(t-1)….e(t-5) where e(i) is the difference between
            # the moving average at ith instant and actual value.
            # q – The lag value where the ACF chart crosses the upper confidence
            # interval for the first time.
            uci = 1.96 / np.sqrt(len(y))
            q = np.where(lag_acf <= uci)[0]
            if q.size > 0:
                q = q[0] + 1
            else:
                q = 1

            # Number of Differences (d): These are the number of nonseasonal
            # differences, i.e. in this case we took the first order difference. So
            # either we can pass that variable and put d=0 or pass the original
            # variable and put d=1. Both will generate same results.
            d = 0

            # Fit model
            model = ARIMA(y, order=(1, 0, 1))
            model_fit = model.fit(disp=False)
            # Make prediction and calculate residual error
            yhat = model_fit.predict(1, len(y))
            prediction_error = y - yhat

            # Write data to file
            pd.DataFrame(prediction_error).to_csv(R_input_filename, index=None,
                                                  header=None)
            # Perform rest of analysis in R
            cmd = '"%s" --no-save "%s" "%s" > "%s" 2>&1' % (
            PARAMETERS['R_path'], R_script, R_input_filename, R_log)
            exit_code = subprocess.call(cmd, shell=True)
            if exit_code != 0:
                raise Exception('Encountered an error while running R')
            # Read results
            res = pd.read_csv(R_input_filename, header=None).values.ravel()
            if len(res) == 1 and res[0] == 'None':
                res = []
            else:
                res = res.tolist()
        except ValueError:
            res = []

    # Store results
    RESULTS['class_comparisons']['change_points'][
        cl] = res  # Pass 1-indexing indices used in R to Matlab

# Show change-points, if user requested to calculate them
if PARAMETERS['classification_change_point_detection']:
    for class_idx, cl in enumerate(class_names):
        if cl == 'SP':
            continue
        # Get timestamps
        res = np.array(RESULTS['class_comparisons']['change_points'][cl]) - 1
        if res.size == 0:
            continue
        # Mark these timestamps
        for i in range(len(res)):
            ax[0, class_idx].axvline(res[i], color='r', lw=2, linestyle='--')
            ax[1, class_idx].axvline(res[i], color='r', lw=2, linestyle='--')
            ax[2, class_idx].axvline(res[i], color='r', lw=2, linestyle='--')
