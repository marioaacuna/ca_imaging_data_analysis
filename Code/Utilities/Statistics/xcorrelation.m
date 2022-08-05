function CCG = xcorrelation(ts1, ts2, lags)

CCG = NaN(length(lags), 1);
for ibin = 1:length(lags)
    CCG(ibin) = corr(ts1, circshift(ts2, lags(ibin)), 'Type','Pearson');
end
