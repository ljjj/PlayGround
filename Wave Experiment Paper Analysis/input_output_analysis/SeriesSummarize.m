function Series = SeriesSummarize(Series)
% This function finishes up final analysis of Series after all data is appended.
% Series : struct defined in SeriesInitialization.m

% order statistics
Series.mean     = nanmean (Series.timeseries,    2); % save the mean
Series.std      = nanstd  (Series.timeseries, 0, 2); % save the std
Series.skewness = skewness(Series.timeseries, 1, 2); % save the skewness
Series.kurtosis = kurtosis(Series.timeseries, 1, 2); % save the kurtosis

% perform PCA
% PCA makes sense because each observation is a PDF and any linear
% combination of sampled points on the PDF represents a general linear
% filter on the PDF, and thus can be any statistics of the random variable
[Series.pca_results.coeff, Series.pca_results.score, ...
    Series.pca_results.latent, Series.pca_results.tsquared, ...
    Series.pca_results.explained, Series.pca_results.mu] = pca(Series.n);
end