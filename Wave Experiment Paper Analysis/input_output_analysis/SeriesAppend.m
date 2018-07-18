function Series = SeriesAppend(Series, timeseries)
% This function appends timeseries to Series.
% Series     : struct defined in SeriesInitialization.m
% timeseries : array 1 x N_times of time Series to append

Series.timeseries = [Series.timeseries; timeseries(:)'];
psd = periodogram(timeseries,rectwin(length(timeseries)),length(timeseries),1/Series.dt);
Series.psd = [Series.psd; psd(:)'];
n = hist(timeseries, Series.ctrs);
n = n/nansum(n)/Series.delta; % normalize to pdf
Series.n = [Series.n; n(:)']; % save the pdf
acf = autocorr(timeseries, Series.N_acf); % find autocorrelation
Series.acf = [Series.acf; acf(:)'];
pacf = parcorr(timeseries, Series.N_acf); % find partial autocorrelation
Series.pacf = [Series.pacf; pacf(:)'];
end