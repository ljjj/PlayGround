% This file takes the final fits from wave_profile_analysis_final_test.m
% and produce figures to be used in the paper.
% reformat for science

%% load data
load('C:\Users\sh683\Dropbox (Personal)\paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\fig2_rawdata.mat')
load('C:\Users\sh683\Dropbox (Personal)\paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\jun_analysis\final\Gaussian_mixture_fits.mat')
load('C:\Users\sh683\Dropbox (Personal)\paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\fig2_material_May2017\jun_analysis\final\wave_speeds_and_cell_counts.mat')
Nexp = numel(data2);
lr_indices = [20, 30, 40, 50, 60, 80, 100, 150, 200]; % indices shown for linear regressions
lr_colors = {'r','m','b',[0.1,0.7,0.8],'g',[1,0.4,0.1],[0.5,0.3,0.2],[0.6,0.6,0.6],'k'}; % colors shown for linear regressions
n_hist = 100; % number of bins in the histogram
gc_colors = {'r','y','g'}; % colors used for Gaussian components
%% show plots
figureS5 = CreateFigSpace('Letter');
[figSize, panelcord] = GetPanelCord(3,3,figureS5);
Labels = [];
plot_i = 0
fignumx = 1;
fignumy = 1;
for i = [9,6,3]
    % initialization
    cp = data2(i).CellDensity.CellPos;
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
    count_max = round(total_count*1.6);
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
    
    % plot linear regressions
    plot_i = plot_i + 1;
%     subplot(3,3,plot_i)
    fignumx = 1;
    hS5A = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
%     hS5A.YAxis.Exponent = 4;
    hold on
    for k = 1:numel(lr_indices) % choose a few places in the wave to show linear regressions
        t0 = th_time(:,lr_indices(k));
        x0 = th_pos(:,lr_indices(k));
        condition = t0 >= 0; % show all
        if i == 1 || i == 3 || i == 18 % wave slowed down after 850 s in these experiments
            condition = t0 < 850;
        end
        % remove data outside bound
        t = t0(condition);
        x = x0(condition);
        n = numel(t);
        
        if n > 5 % fit only with enough points
            % least squares fitting
            xt = mean(x.*t);
            tt = mean(t.^2);
            b1 = (xt - mean(x)*mean(t))/(tt - mean(t)^2); % slope
            b0 = mean(x) - b1*mean(t); % intercept
            plot(t0/60, x0/1000, 'o', 'color', lr_colors{k}) % show data
            plot(t0/60, (b1*t0 + b0)/1000, '--', 'color', lr_colors{k}) % show fit
        end
    end
    Labels(end+1) = xlabel('Time (min)')
    Labels(end+1) = ylabel('Position (mm)')
    axis tight
    ylim([0, 10])
    
    % plot regressed speed vs. cumulative counts
    plot_i = plot_i + 1;
%     subplot(3,3,plot_i)
    fignumx = 2;
    hS5B = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
%     hS5B.XAxis.Exponent = 4;
    hold on
    plot(count_thresh{i}./10000, (speeds_ls{i}).*60./1000, 'k') % mean speed
    plot(count_thresh{i}./10000, (speeds_ls{i}+se_speeds_ls{i}).*60./1000, 'k:')
    plot(count_thresh{i}./10000, (speeds_ls{i}-se_speeds_ls{i}).*60./1000, 'k:')
    for k = 1:numel(lr_indices) % show where along the wave the regressions are plotted
        scatter((count_thresh{i}(lr_indices(k)))./10000, (speeds_ls{i}(lr_indices(k))).*60./1000, [], lr_colors{k}, 'o', 'filled')
    end
    plot(total_cells(i)./10000 * [1 1], [-1 7], 'k')
    plot(total_cells_l(i)./10000 * [1 1], [-1 7], 'k--')
    plot(total_cells_h(i)./10000 * [1 1], [-1 7], 'k--')
    Labels(end+1) = xlabel('n''th cell from the front (10^4)')
    Labels(end+1) = ylabel('Speed (mm/min)')
    axis tight
    ylim([0, 0.4])

    % plot regressed speed histogram and Gaussian mixture fit
    xs = -1:0.001:7;
    plot_i = plot_i + 1;
%     subplot(3,3,plot_i)
    fignumx = 3;
    hS5C = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
    h = histogram((speeds_ls{i}).*60./1000, n_hist, 'Normalization', 'pdf', 'Orientation', 'horizontal');
    hold on
    for j = 1:N_G_choices(1) % loop over Gaussian components
        plot((props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2)).*1000/60,...
            (xs).*60./1000, gc_colors{j}, 'LineWidth', 1)
    end
    Labels(end+1) = xlabel('pdf')
    Labels(end+1) = ylabel('Speed (mm/min)')
    xlim([0, 40])
    ylim([0, 0.4])
    fignumy = fignumy+1;
end
prettyFig('lw',1.5,'fs',10,'x_minor_ticks','true','y_minor_ticks','true','tick_dir','in','axis_box','off','tick_length',0.025)
for ll = 1:length(Labels)
    set(Labels(ll),'FontSize',12)
end
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','A','units','centimeters','Position',[-1  figSize.y+14.25 0]);
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','delete \mu','units','centimeters','Position',[-1  figSize.y+15.25 0]);

text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','B','units','centimeters','Position',[-1  figSize.y+14.25 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','C','units','centimeters','Position',[-1  figSize.y+14.25 0]);
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','D','units','centimeters','Position',[-1  figSize.y+7.25 0]);
text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','E','units','centimeters','Position',[-1  figSize.y+7.25 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','F','units','centimeters','Position',[-1  figSize.y+7.25 0]);
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','G','units','centimeters','Position',[-1  figSize.y+0.25 0]);
text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','H','units','centimeters','Position',[-1  figSize.y+0.25 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','I','units','centimeters','Position',[-1  figSize.y+0.25 0]);