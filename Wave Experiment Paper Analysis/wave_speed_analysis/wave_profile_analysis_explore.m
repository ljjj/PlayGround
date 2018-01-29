% This file explores how to fit linear regressions and Gaussian mixture
% models, trying out many possibilities. It also generates some exploratory
% plots.
%% load data
load('C:\Users\jl2345\Dropbox (emonetlab)\paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\fig2_rawdata.mat')
%% find wave speed with threshold cell counts
plot_all_counts_profile = 0;
plot_profiles_sweep = 0;
plot_all_wavefront_series = 1;
delta_x = 122; % uncertainty in space due to image bin size when counting cells, um
delta_t = 3.2; % uncertainty in time due to image taking frequency, s
delta = (delta_x / delta_t)^2; % delta >> 1 makes Deming regression = OLS
Nexp = numel(data2);
count_thresh = cell(1, Nexp);
speeds_ls = cell(1, Nexp);
se_speeds_ls = cell(1, Nexp);
speeds_dm = cell(1, Nexp);
for i = 1:Nexp % different experiments
    % initialization
    if isempty(data2(i).CellDensity)
        continue
    end
    t_str = ['Asp = ',num2str(data2(i).AspConc),' \muM'];
    if isempty(data2(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    cp = data2(i).CellDensity.CellPos;
    if isempty(cp(2).pos)
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    N_s = numel(sweeps);
    
    % test plot all profiles
    if plot_all_counts_profile
        figure('visible','off')
        hold on
        xlabel('x (\mum)')
        ylabel('t (s)')
        zlabel('cell counts')
    end
    
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
    count_max = round(total_count*1.6);
    dc = round(count_max/300);
    count_thresh{i} = dc:dc:count_max; % use cell count threshold to determine wave front
    N_t = numel(count_thresh{i});
    th_pos = nan(N_s, N_t);
    th_time = nan(N_s, N_t);
    
    for j = 1:N_s % different sweeps
        if plot_all_counts_profile
            plot3(cp(sweeps(j)).pos, cp(sweeps(j)).time, cp(sweeps(j)).counts)
        end
        % test plot single time profile
        if plot_profiles_sweep
            figure('visible','off')
            subplot(2,1,1)
            plot(cp(sweeps(j)).pos,cp(sweeps(j)).counts)
            xlabel('x (\mum)')
            ylabel('cell counts')
            title([t_str, ' t = ', num2str(cp(sweeps(j)).time(1)/60), ' min'])
            subplot(2,1,2)
            plot(cp(sweeps(j)).pos,cumsum(cp(sweeps(j)).counts,'reverse'))
            xlabel('x (\mum)')
            ylabel('cumulative cell counts')
            saveas(gcf, ['profiles_',num2str(i),'_sweep_',num2str(j),'.fig'])
            saveas(gcf, ['profiles_',num2str(i),'_sweep_',num2str(j),'.jpg'])
            close(gcf)
        end
        
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
    
    if plot_all_counts_profile
        view(20,30)
        title(t_str)
        saveas(gcf, ['all_counts_profiles_',num2str(i),'.fig'])
        saveas(gcf, ['all_counts_profiles_',num2str(i),'.jpg'])
        close(gcf)
    end
    
    % find wave speed and its se
%     x_bound = [0, 8000]; % only fit within this interval
    speeds_ls{i} = nan(1,N_t);
    se_speeds_ls{i} = nan(1,N_t);
    speeds_dm{i} = nan(1,N_t);
    % test plot wave front position-time
    if plot_all_wavefront_series
        figure('visible','off')
        hold on
    end
    for k = 1:N_t
        t0 = th_time(:,k);
        x0 = th_pos(:,k);
        condition = t0 >= 0; % all
        if i == 1 || i == 3 || i == 18
            condition = t0 < 850;
        end
%         condition = x0 > x_bound(1) & x0 < x_bound(2);
        if plot_all_wavefront_series && mod(k,20)==0
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
            if plot_all_wavefront_series && mod(k,20)==0
                plot(t0, b1*t0 + b0, '--')
            end
            e2 = sum((x - b1*t - b0).^2);
            se_speeds_ls{i}(k) = sqrt(e2/(n-2)/n)/std(x);

            % Deming regression
            xm = mean(x);
            tm = mean(t);
            sxx = sum((x-xm).^2)/(n-1);
            sxt = sum((x-xm).*(t-tm))/(n-1);
            stt = sum((t-tm).^2)/(n-1);
            b1 = (sxx - delta*stt + sqrt((sxx - delta*stt)^2 + 4*delta*sxt^2)) / (2*sxt);
            b0 = mean(x) - b1*mean(t);
            speeds_dm{i}(k) = b1;
            if plot_all_wavefront_series && mod(k,20)==0
                plot(t0, b1*t0 + b0)
            end
        end
    end
    
    if plot_all_wavefront_series
        xlabel('t (s)')
        ylabel('x (\mum)')
        title(t_str)
        saveas(gcf, ['all_wavefront_series_',num2str(i),'.fig'])
        saveas(gcf, ['all_wavefront_series_',num2str(i),'.jpg'])
        close(gcf)
    end
end
%% initialize Gaussian mixture fits
N_G_choices = [2, 3];
mus = zeros(Nexp, sum(N_G_choices));
sigmas = zeros(Nexp, sum(N_G_choices));
props = zeros(Nexp, sum(N_G_choices));
%% fit Gaussian mixture model to the histogram of speeds
plot_all_threshold_regression_fits = 1;
n_hist = 50; % number of bins in the histogram
colors = {'r', 'g', 'r', 'y', 'g'}; % colors used for Gaussian components
for i = 18 % sometimes convergence fails and needs to try again
    if isempty(speeds_ls{i})
        continue
    end
    t_str = ['Asp = ',num2str(data2(i).AspConc),' \muM'];
    if isempty(data2(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    
    % fit Gaussian mixture model of N_G
    j = 0; % saving index
    for N_G = N_G_choices
        options = statset('MaxIter',10000);
        X = speeds_ls{i}';
        X(isnan(X)) = [];
        X(X<0) = [];
%         y = 2*ones(size(X,1),1);
%         y(X > max(X)-1) = N_G;
%         y(X < 1) = 1;
        GMModel = fitgmdist(X, N_G, 'Options', options);
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
    
    % test plot regression results and fits of histograms
    if plot_all_threshold_regression_fits
%         figure('visible','off')
        % plot regressed speed vs. cumulative counts
        subplot(1,3,1)
        hold on
        plot(count_thresh{i}, speeds_ls{i}, 'k')
        plot(count_thresh{i}, speeds_ls{i}+se_speeds_ls{i}, 'k--')
        plot(count_thresh{i}, speeds_ls{i}-se_speeds_ls{i}, 'k--')
        xlabel('cumulative cell #')
        ylabel('speed (\mum/s)')
        axis tight
        ylim([-1,7])
        title(t_str)
        
        % plot regressed speed histogram and Gaussian mixture fit
        xs = -1:0.01:7;
        subplot(1,3,2)
        hold off
        h = histogram(speeds_ls{i}, n_hist, 'Normalization', 'pdf', 'Orientation', 'horizontal');
        hold on
        for j = 1:N_G_choices(1)
            plot(props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2),...
                xs, colors{j}, 'LineWidth', 1)
        end
        xlabel('pdf')
        ylabel('speed (\mum/s)')
        title('2 Gaussian Mixture Fit')
        xlim([0,2])
        ylim([-1,7])
        
        subplot(1,3,3)
        hold off
        h = histogram(speeds_ls{i}, n_hist, 'Normalization', 'pdf', 'Orientation', 'horizontal');
        hold on
        for j = (N_G_choices(1)+1):(N_G_choices(1)+N_G_choices(2))
            plot(props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2),...
                xs, colors{j}, 'LineWidth', 1)
        end
        xlabel('pdf')
        ylabel('speed (\mum/s)')
        title('3 Gaussian Mixture Fit')
        xlim([0,2])
        ylim([-1,7])
        
        % save plot
%         set(gcf, 'PaperPosition', [0 0 13 4]);
%         saveas(gcf, ['all_threshold_regression_fits_',num2str(i),'.fig'])
%         saveas(gcf, ['all_threshold_regression_fits_',num2str(i),'.jpg'])
%         close(gcf)
    end
end
%% save Gaussian mixture fits, don't run this without know what you did
% save('Gaussian_mixture_fits.mat','mus','sigmas','props','N_G_choices','speeds_ls','se_speeds_ls','speeds_dm')
%% load Gaussian mixture fits
load('Gaussian_mixture_fits.mat')
%% finalize wave speed
speed_ind = 5*ones(1, Nexp);
% speed_ind([13 16 19]) = 2;
wave_speed = nan(1, Nexp);
wave_speed_se = nan(1, Nexp);
for i = 1:Nexp
    if ~isempty(speeds_ls{i})
        wave_speed(i) = mus(i, speed_ind(i));
        wave_speed_se(i) = sigmas(i, speed_ind(i));
    end
end
%% find where the wave ends
plot_wave_number_selection = 1;
total_cells = nan(1, Nexp);
for i = 1:Nexp
    if isempty(speeds_ls{i})
        continue
    end
    t_str = ['Asp = ',num2str(data2(i).AspConc),' \muM'];
    if isempty(data2(i).TBdist)
        t_str = [t_str, ' no P(TB)'];
    else
        t_str = [t_str, ' has P(TB)'];
    end
    
    if speed_ind(i) == 5
%         wave_end_speed_l = max(mus(i,4) - sigmas(i,4), mus(i,3) + 3*sigmas(i,3));
%         wave_end_speed_h = mus(i,4);
        ii = 4;
        si = 0;
    else
%         wave_end_speed_l = mus(i,1);
%         wave_end_speed_h = mus(i,1) + sigmas(i,1);
        ii = 1;
        si = 0;
    end
%     wave_end_range = speeds_ls{i} > wave_end_speed_l & speeds_ls{i} < wave_end_speed_h;
    wave_end_range = speeds_ls{i} < mus(i,ii) + sigmas(i,ii) * si; % cut off at the back where the speed is in the ii component, with si*std away from mean
    [~, imax] = max(speeds_ls{i});
    wave_end_range(1:imax) = 0; % not cut off in front of the wave
    for j = 1:numel(speeds_ls{i})
        if wave_end_range(j)
            total_cells(i) = count_thresh{i}(j);
            break
        end
    end
    if plot_wave_number_selection
        figure('visible','off')
        hold on
        plot(count_thresh{i}, speeds_ls{i}, 'k')
        plot(count_thresh{i}, speeds_ls{i}+se_speeds_ls{i}, 'k--')
        plot(count_thresh{i}, speeds_ls{i}-se_speeds_ls{i}, 'k--')
        plot(total_cells(i) * [1 1], [-1 7], 'r-.')
        xlabel('cumulative cell #')
        ylabel('speed (\mum/s)')
        ylim([-1,7])
        title(t_str)
        saveas(gcf, ['wave_number_selection_',num2str(i),'.fig'])
        saveas(gcf, ['wave_number_selection_',num2str(i),'.jpg'])
        close(gcf)
    end
end
%% save final results
save('wave_speeds_and_cell_counts.mat','wave_speed','wave_speed_se','total_cells')