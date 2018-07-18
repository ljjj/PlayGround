function mi = SeriesMutualInformation(Series1, Series2, lag)
% This function calculates the mutual information between the last time
% series of two Series structs in bits.
% Series1,2  : structs defined in SeriesInitialization.m
% lag        : lag between the time series
%              lag < 0 means Series1 causes Series2

% get the two time series
ts1 = Series1.timeseries(end,:);
ts2 = Series2.timeseries(end,:);

% calculate joint distribution
ts2_lagged = lagmatrix(ts2, lag);
ts = [ts1(:) ts2_lagged(:)];
p12 = hist3(ts, {Series1.ctrs, Series2.ctrs});
p12 = p12/nansum(nansum(p12))/Series1.delta/Series2.delta;

% calculate individual distributions
% p1 = hist(ts1, Series1.ctrs);
% p1 = p1/nansum(p1)/Series1.delta; % normalize to pdf
% p2 = hist(ts2, Series2.ctrs);
% p2 = p2/nansum(p2)/Series2.delta; % normalize to pdf
p1 = sum(p12,2)*Series2.delta;
p2 = sum(p12,1)*Series1.delta;

% check consistency
% if any(sum(p12,2)*Series2.delta ~= p1(:)) || any(sum(p12,1)*Series1.delta ~= p2(:)')
%     error('distribution calculation incorrect')
% end

% calculate mutual informatoin
p1 = repmat(p1, 1, Series2.N_bins);
p2 = repmat(p2, Series1.N_bins, 1);
mi = nansum(nansum(p12.*log(p12./p1./p2)))*Series1.delta*Series2.delta/log(2);
end