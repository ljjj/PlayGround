% This file builds on wave_profile_analysis_explore.m, removing redundant
% plots and focus on finalizing Gaussian mixture routine by checking the
% regularization value. This file generates the final fits.
%% load data
load('C:\Users\jl2345\Dropbox (emonetlab)\users\setsu_kato\papers\xiongfei_setsu\manuscript\Figures\codes for figs\raw data\data_aTc20_400.mat')
%% initialization
Nexp = numel(data_aTc20_400);
% linear fit
delta_x = 130; % uncertainty in space due to image bin size when counting cells, um
delta_t = 1.5; % uncertainty in time due to image taking frequency, s
delta = (delta_x / delta_t)^2; % delta >> 1 makes Deming regression = OLS
count_thresh = cell(1, Nexp);
speeds_ls = cell(1, Nexp);
se_speeds_ls = cell(1, Nexp);

% Gaussian mixture
N_G_choices = 3;
mus = zeros(Nexp, sum(N_G_choices));
sigmas = zeros(Nexp, sum(N_G_choices));
props = zeros(Nexp, sum(N_G_choices));
n_hist = 100; % number of bins in the histogram
colors = {'r', 'y', 'g'}; % colors used for Gaussian components
reg = 0.01;
%% fits
for i = 1:Nexp % different experiments
    % initialization
    if isempty(data_aTc20_400(i).CellDensity)
        continue
    end
    t_str = ['Asp = ',num2str(data_aTc20_400(i).AspConc),' \muM'];
    if isempty(data_aTc20_400(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    cp = data_aTc20_400(i).CellDensity.CellPos;
    if isempty(cp(2).pos)
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    N_s = numel(sweeps);
    
    % find approximate number of cells in the wave
    cm = cp(5).counts_smooth; % later wave for accuracy but not too late to keep more cells
    dsd = diff(sign( diff(cm) ));
    den_min_ind = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
    if isempty(den_min_ind)
        den_min_ind = 1;
    end
    den_min_ind = den_min_ind(1);
    cc = cumsum(cp(5).counts, 'reverse');
    total_count = cc(den_min_ind);
    % used to determine the maximal threshold to use
    count_max = min(round(total_count*1.7),max(cc)*0.9);
    dc = round(count_max/300);
    count_thresh{i} = dc:dc:count_max; % use cell count threshold to determine wave front
    N_t = numel(count_thresh{i});
    th_pos = nan(N_s, N_t);
    th_time = nan(N_s, N_t);
    
    for j = 1:N_s % different sweeps
        cum_counts = cumsum(cp(sweeps(j)).counts, 'reverse');
        for k = 1:N_t % different wave front threshold choices
            ind = find(cum_counts > count_thresh{i}(k));
            if isempty(ind)
                l = 1;
            else
                l = ind(end); % index where cell count exceeds threshold
            end
            th_pos(j,k) = cp(sweeps(j)).pos(l);
            th_time(j,k) = cp(sweeps(j)).time(l);
        end
    end
    
    % find wave speed and its se
    speeds_ls{i} = nan(1,N_t);
    se_speeds_ls{i} = nan(1,N_t);
    figure('visible','off')
    subplot(1,3,1)
    hold on
    for k = 1:N_t
        t0 = th_time(:,k);
        x0 = th_pos(:,k);
        condition = t0 >= 0; % all
        if i == 2 || i == 4
            condition = t0 < 1000;
        end
        if mod(k,20)==0
            plot(t0, x0, 'o')
        end
        % remove data outside bound
        t = t0(condition);
        x = x0(condition);
        n = numel(t);
        
        if n > 5 % fit only with enough points
            % least squares fitting
            xt = mean(x.*t);
            tt = mean(t.^2);
            b1 = (xt - mean(x)*mean(t))/(tt - mean(t)^2);
            b0 = mean(x) - b1*mean(t);
            speeds_ls{i}(k) = b1;
            if mod(k,20)==0
                plot(t0, b1*t0 + b0, '--')
            end
            e2 = sum((x - b1*t - b0).^2);
            se_speeds_ls{i}(k) = sqrt(e2/(n-2)/n)/std(x);
        end
    end
    xlabel('t (s)')
    ylabel('x (\mum)')
    
    % plot regressed speed vs. cumulative counts
    subplot(1,3,2)
    hold on
    plot(count_thresh{i}, speeds_ls{i}, 'k')
    plot(count_thresh{i}, speeds_ls{i}+se_speeds_ls{i}, 'k--')
    plot(count_thresh{i}, speeds_ls{i}-se_speeds_ls{i}, 'k--')
    xlabel('cumulative cell #')
    ylabel('speed (\mum/s)')
    axis tight
    ylim([0,7])
    title(t_str)

    % fit Gaussian mixture model of N_G
    j = 0; % saving index
    for N_G = N_G_choices
        options = statset('MaxIter',10000,'TolFun',1e-10);
        X = speeds_ls{i}';
        X(isnan(X)) = [];
        X(X<0) = [];
        GMModel = fitgmdist(X, N_G, 'RegularizationValue', reg, 'Replicates', 10, 'Options', options);
        mu = GMModel.mu;
        sigma = squeeze(GMModel.Sigma);
        prop = GMModel.ComponentProportion;
        % order the models
        [mu, sort_ind] = sort(mu);
        sigma = sigma(sort_ind);
        prop = prop(sort_ind);
        % save
        mus(i, (j+1):(j+N_G)) = mu(:);
        sigmas(i, (j+1):(j+N_G)) = sigma(:);
        props(i, (j+1):(j+N_G)) = prop(:);
        j = j + N_G; % increase the counter offset
    end
    
    % plot regressed speed histogram and Gaussian mixture fit
    xs = -1:0.001:7;
    subplot(1,3,3)
    hold off
    h = histogram(speeds_ls{i}, n_hist, 'Normalization', 'pdf', 'Orientation', 'horizontal');
    hold on
    for j = 1:N_G_choices(1)
        plot(props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2),...
            xs, colors{j}, 'LineWidth', 1)
    end
    xlabel('pdf')
    ylabel('speed (\mum/s)')
    xlim([0,2])
    ylim([0,7])

    % save plot
    set(gcf, 'PaperPosition', [0 0 15 4]);
    saveas(gcf, ['all_threshold_regression_fits_',num2str(i),'_reg_',num2str(reg),'.jpg'])
    close(gcf)
end
%% save Gaussian mixture fits, don't run this without know what you did
save('..\Gaussian_mixture_fits.mat','mus','sigmas','props','N_G_choices','speeds_ls','se_speeds_ls','reg')
%% load Gaussian mixture fits
load('..\Gaussian_mixture_fits.mat')
%% finalize wave speed
wave_speed = nan(1, Nexp);
wave_speed_sd = nan(1, Nexp);
for i = 1:Nexp
    if ~isempty(speeds_ls{i})
        wave_speed(i) = mus(i, 3);
        wave_speed_sd(i) = sigmas(i, 3);
    end
end
%% find where the wave ends
total_cells = nan(1, Nexp);
total_cells_l = nan(1, Nexp);
total_cells_h = nan(1, Nexp);
for i = 1:Nexp
    if isempty(speeds_ls{i})
        continue
    end
    t_str = ['Asp = ',num2str(data_aTc20_400(i).AspConc),' \muM'];
    if isempty(data_aTc20_400(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    
    % index to define where the wave is
    wave_end = speeds_ls{i} < mus(i,2);
    wave_end_l = speeds_ls{i} < mus(i,2) - sigmas(i,2);
    wave_end_h = speeds_ls{i} < mus(i,2) + sigmas(i,2);
    [~, imax] = max(speeds_ls{i});
    
    % not in the front of the wave
    wave_end(1:imax) = 0;
    wave_end_l(1:imax) = 0;
    wave_end_h(1:imax) = 0;
    
    for j = 1:numel(speeds_ls{i})
        if wave_end(j)
            total_cells(i) = count_thresh{i}(j);
            break
        end
    end
    
    for j = 1:numel(speeds_ls{i})
        if wave_end_l(j)
            total_cells_l(i) = count_thresh{i}(j);
            break
        end
    end
    
    for j = 1:numel(speeds_ls{i})
        if wave_end_h(j)
            total_cells_h(i) = count_thresh{i}(j);
            break
        end
    end
    
    figure('visible','off')
    hold on
    plot(count_thresh{i}, speeds_ls{i}, 'k')
    plot(count_thresh{i}, speeds_ls{i}+se_speeds_ls{i}, 'k--')
    plot(count_thresh{i}, speeds_ls{i}-se_speeds_ls{i}, 'k--')
    plot(total_cells(i) * [1 1], [-1 7], 'r-')
    plot(total_cells_l(i) * [1 1], [-1 7], 'r-.')
    plot(total_cells_h(i) * [1 1], [-1 7], 'r-.')
    xlabel('cumulative cell #')
    ylabel('speed (\mum/s)')
    ylim([-1,7])
    title(t_str)
    saveas(gcf, ['wave_number_selection_',num2str(i),'.fig'])
    saveas(gcf, ['wave_number_selection_',num2str(i),'.jpg'])
    close(gcf)
end
%% save final results
save('..\wave_speeds_and_cell_counts.mat','wave_speed','wave_speed_sd','total_cells','total_cells_l','total_cells_h')