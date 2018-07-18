function Series = SeriesInitialization(name, N_bins, min, max, dt, N_times, N_acf)
% This function initializes Series struct.
% Series : struct with fields:
%          name       : Series name
%          N_bins     : number of bins in histogram
%          min        : histogram lower bound
%          max        : histogram upper bound
%          delta      : histogram bin size
%          ctrs       : array 1 x N_bins of histogram centers
%          timeseries : array N_tracks x N_times of time Series
%          time       : array 1 x N_times of time
%          dt         : time interval
%          freq       : array 1 x floor(N_times/2)+1 of frequency
%          psd        : array 1 x floor(N_times/2)+1 of power spectral density
%          n          : array N_tracks x N_bins of histogram pdfs
%          mean       : array N_tracks x 1 of means
%          std        : array N_tracks x 1 of stds
%          skewness   : array N_tracks x 1 of skewnesses
%          kurtosis   : array N_tracks x 1 of kurtoses
%          N_acf      : autocorrelation maximal lag
%          acf        : autocorrelation values
%          pacf       : partial-autocorrelation values
%          acf_lags   : array 1 x N_acf+1 of lag times for acf and pacf
%          pca_results: struct of outputs from pca

Series = [];
Series.name = name;
Series.N_bins = N_bins;
Series.min = min;
Series.max = max;
Series.delta = (Series.max - Series.min)/Series.N_bins;
Series.ctrs = (Series.min+Series.delta/2):Series.delta:(Series.max-Series.delta/2);
Series.timeseries = [];
Series.time = (0:(N_times-1))*dt;
Series.dt = dt;
Series.freq = 0:1/(N_times*dt):1/(2*dt);
Series.psd = [];
Series.n = [];
Series.mean = [];
Series.std = [];
Series.skewness = [];
Series.kurtosis = [];
Series.N_acf = N_acf;
Series.acf = [];
Series.pacf = [];
Series.acf_lags = (0:N_acf)*dt;
Series.pca_results = [];
end