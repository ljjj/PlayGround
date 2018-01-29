% This file takes the final fits from wave_profile_analysis_final_test.m
% and produce figures to be used in the paper.
%% load data
load('C:\Users\jl2345\Dropbox (emonetlab)\users\setsu_kato\papers\xiongfei_setsu\manuscript\Figures\codes for figs\raw data\data_aTc20_400.mat')
load('Gaussian_mixture_fits.mat')
load('wave_speeds_and_cell_counts.mat')
Nexp = numel(data_aTc20_400);
lr_indices = [20, 30, 40, 50, 60, 80, 100, 150, 200]; % indicies shown for linear regressions
lr_colors = {'r','m','b',[0.1,0.7,0.8],'g',[1,0.4,0.1],[0.5,0.3,0.2],[0.6,0.6,0.6],'k'}; % colors shown for linear regressions
n_hist = 100; % number of bins in the histogram
gc_colors = {'r','y','g'}; % colors used for Gaussian components
count_thresh = cell(1, Nexp);
%% show plots
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
    
    figure('visible','off')
    % plot linear regressions
    subplot(1,3,1)
    hold on
    for k = 1:numel(lr_indices)
        t0 = th_time(:,lr_indices(k));
        x0 = th_pos(:,lr_indices(k));
        condition = t0 >= 0; % all
        if i == 2 || i == 4
            condition = t0 < 1000;
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
            plot(t0, x0, 'o', 'color', lr_colors{k})
            plot(t0, b1*t0 + b0, '--', 'color', lr_colors{k})
        end
    end
    xlabel('t (s)')
    ylabel('x (\mum)')
    axis tight
    ylim([0, 10000])
    
    % plot regressed speed vs. cumulative counts
    subplot(1,3,2)
    hold on
    plot(count_thresh{i}, speeds_ls{i}, 'k')
    plot(count_thresh{i}, speeds_ls{i}+se_speeds_ls{i}, 'k:')
    plot(count_thresh{i}, speeds_ls{i}-se_speeds_ls{i}, 'k:')
    for k = 1:numel(lr_indices)
        scatter(count_thresh{i}(lr_indices(k)), speeds_ls{i}(lr_indices(k)), [], lr_colors{k}, 'o', 'filled')
    end
    plot(total_cells(i) * [1 1], [-1 7], 'k')
    plot(total_cells_l(i) * [1 1], [-1 7], 'k--')
    plot(total_cells_h(i) * [1 1], [-1 7], 'k--')
    xlabel('cumulative cell #')
    ylabel('speed (\mum/s)')
    axis tight
    ylim([0, 7])
    title(t_str)

    % plot regressed speed histogram and Gaussian mixture fit
    xs = -1:0.001:7;
    subplot(1,3,3)
    hold off
    h = histogram(speeds_ls{i}, n_hist, 'Normalization', 'pdf', 'Orientation', 'horizontal');
    hold on
    for j = 1:N_G_choices(1)
        plot(props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2),...
            xs, gc_colors{j}, 'LineWidth', 1)
    end
    xlabel('pdf')
    ylabel('speed (\mum/s)')
    xlim([0, 2])
    ylim([0, 7])

    % save plot
    set(gcf, 'PaperPosition', [0 0 15 4]);
    print(gcf, ['exp_',num2str(i),'.jpg'], '-djpeg', '-r300', '-opengl')
    close(gcf)
end