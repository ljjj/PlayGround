% This file takes the final fits from wave_profile_analysis_final_test.m
% and produce figures to be used in the paper.
% reformat for science

%% load data
load('..\..\fig2_rawdata.mat')
load('..\..\wave_speed_analysis\final\Gaussian_mixture_fits.mat')
load('..\..\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat')
load('..\..\align_density_profile_analysis\explorations\5\4\results.mat')
Nexp = numel(data2);
lr_indices = [20, 30, 40, 50, 60, 80, 100, 150, 200]; % indices shown for linear regressions
lr_colors = {'r','m','b',[0.1,0.7,0.8],'g',[1,0.4,0.1],[0.5,0.3,0.2],[0.6,0.6,0.6],'k'}; % colors shown for linear regressions
n_hist = 100; % number of bins in the histogram
gc_colors = {'r','y','g'}; % colors used for Gaussian components
%% show plots
figureS5 = CreateFigSpace('Letter');
[figSize, panelcord] = GetPanelCord(5,3,figureS5);
panelcord.y = (panelcord.y - 2) * 0.84 + 2;
figSize.y = figSize.y * 0.8;
Labels = [];
fignumx = 1;
fignumy = 1;
for i = [7 17 3]
    % initialization
    cp = data2(i).CellDensity.CellPos;
    if isempty(cp(2).pos)
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    nSweeps = numel(sweeps);
    sweepTimes = cellfun(@mean, {cp(sweeps).time});
    
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
    nThresholds = numel(count_thresh{i});
    th_pos = nan(nSweeps, nThresholds);
    th_time = nan(nSweeps, nThresholds);
    
    for j = 1:nSweeps % different sweeps
        cum_counts = cumsum(cp(sweeps(j)).counts, 'reverse');
        for k = 1:nThresholds % different wave front threshold choices
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
    
    % extract KDE results
    sampleZ = vertcat(cp(sweeps).pos) - finalSpeed(i)*vertcat(cp(sweeps).time) - finalShift(i);
    sampleLabels = finalLabels{i};
    sampleCounts = finalCounts{i};
    waveCounts = sampleLabels.*sampleCounts; % cell counts in wave component, sweep weight doesn't matter for each sweep is independent
    
    % plot linear regressions
    fignumy = 1;
    hS5A = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
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
    Labels(end+1) = xlabel('Time (min)');
    Labels(end+1) = ylabel('Position (mm)');
    xlim([0 15])
    ylim([0 10])
    
    % plot regressed speed vs. cumulative counts
    fignumy = 2;
    hS5B = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
    plot(speeds_ls{i}*60/1000, count_thresh{i}/10000, 'k') % mean speed
    plot((speeds_ls{i}+se_speeds_ls{i})*60/1000, count_thresh{i}/10000, 'k:')
    plot((speeds_ls{i}-se_speeds_ls{i})*60/1000, count_thresh{i}/10000, 'k:')
    for k = 1:numel(lr_indices) % show where along the wave the regressions are plotted
        scatter(speeds_ls{i}(lr_indices(k))*60/1000, count_thresh{i}(lr_indices(k))/10000, [], lr_colors{k}, 'o', 'filled')
    end
    plot([0, 0.4], total_cells(i)/10000 * [1 1], 'k')
    plot([0, 0.4], total_cells_l(i)/10000 * [1 1], 'k--')
    plot([0, 0.4], total_cells_h(i)/10000 * [1 1], 'k--')
    Labels(end+1) = xlabel('Speed (mm/min)');
    Labels(end+1) = ylabel('Cumulative cell # (10^4)');
    xlim([0 0.4])
    ylim([0 10])

    % plot regressed speed histogram and Gaussian mixture fit
    xs = -1:0.001:7;
    fignumy = 3;
    hS5C = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
    h = histogram(speeds_ls{i}*60/1000, n_hist, 'Normalization', 'pdf');
    for j = 1:N_G_choices(1) % loop over Gaussian components
        plot(xs*60/1000, props(i,j) * exp(-(xs - mus(i,j)).^2/(2 * sigmas(i,j)^2)) / sqrt(2*pi*sigmas(i,j)^2)*1000/60, gc_colors{j}, 'LineWidth', 1)
    end
    Labels(end+1) = xlabel('Speed (mm/min)');
    Labels(end+1) = ylabel('pdf');
    xlim([0 0.4])
    ylim([0 40])
    
    % plot density profiles at different time points and mean KDE
    fignumy = 4;
    hS5D = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
    for j = 2:nSweeps-5
        sweepColor = j/(nSweeps-5) * [1 0 0] + (nSweeps-5-j)/(nSweeps-5) * [0 0 1];
        rhoSweep = mvksdensity(sampleZ(j,:)',zVec(:),'Bandwidth',bw,'Weights',waveCounts(j,:)');
        plot(zVec/1000, rhoSweep*1e3, 'color', sweepColor)
    end
    plot(zVec/1000, finalRho(:,i)*1e3, 'k--', 'lineWidth', 2)
    Labels(end+1) = xlabel('Moving frame z (mm)');
    Labels(end+1) = ylabel('pdf (10^-^3)');
    xlim([-2 3])
    ylim([0 1])
    
    % plot number of cells in the wave over time
    fignumy = 5;
    hS5E = axes('Units','centimeters','Position',[panelcord.x(fignumx) panelcord.y(fignumy) figSize.x figSize.y]);hold on
    sweepColors = sweepTimes(2:end-5)'/sweepTimes(end-5) * [1 0 0] + (sweepTimes(end-5)-sweepTimes(2:end-5)')/sweepTimes(end-5) * [0 0 1];
    x = sweepTimes(2:end-5)/60;
    y = finalWaveCounts{i}(2:end-5)/10000;
    scatter(x, y, [], sweepColors, 'filled')
    percLoss = (y(1)-y(end))/y(1);
    text(1, 10, [num2str(percLoss*100,2),'% loss over 13 min'])
    Labels(end+1) = xlabel('Time (min)');
    Labels(end+1) = ylabel('Cell # in wave (10^4)');
    xlim([0 15])
    ylim([0 10])
    
    fignumx = fignumx+1;
end
prettyFig('lw',1.5,'plw',1,'fs',10,'x_minor_ticks','true','y_minor_ticks','true','tick_dir','in','axis_box','off','tick_length',0.025)
for ll = 1:length(Labels)
    set(Labels(ll),'FontSize',12)
end
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','a','units','centimeters','Position',[figSize.x - 15.2 figSize.y+0.4 0]);
text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','b','units','centimeters','Position',[figSize.x - 15.2 figSize.y+0.4 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','c','units','centimeters','Position',[figSize.x - 15.2 figSize.y+0.4 0]);
text('Parent',hS5D,'FontWeight','bold','FontSize',16,'String','d','units','centimeters','Position',[figSize.x - 15.2 figSize.y+0.4 0]);
text('Parent',hS5E,'FontWeight','bold','FontSize',16,'String','e','units','centimeters','Position',[figSize.x - 15.2 figSize.y+0.4 0]);
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','f','units','centimeters','Position',[figSize.x - 9.5  figSize.y+0.4 0]);
text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','g','units','centimeters','Position',[figSize.x - 9.5  figSize.y+0.4 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','h','units','centimeters','Position',[figSize.x - 9.5  figSize.y+0.4 0]);
text('Parent',hS5D,'FontWeight','bold','FontSize',16,'String','i','units','centimeters','Position',[figSize.x - 9.5  figSize.y+0.4 0]);
text('Parent',hS5E,'FontWeight','bold','FontSize',16,'String','j','units','centimeters','Position',[figSize.x - 9.5  figSize.y+0.4 0]);
text('Parent',hS5A,'FontWeight','bold','FontSize',16,'String','k','units','centimeters','Position',[figSize.x - 3.8  figSize.y+0.4 0]);
text('Parent',hS5B,'FontWeight','bold','FontSize',16,'String','l','units','centimeters','Position',[figSize.x - 3.8  figSize.y+0.4 0]);
text('Parent',hS5C,'FontWeight','bold','FontSize',16,'String','m','units','centimeters','Position',[figSize.x - 3.8  figSize.y+0.4 0]);
text('Parent',hS5D,'FontWeight','bold','FontSize',16,'String','n','units','centimeters','Position',[figSize.x - 3.8  figSize.y+0.4 0]);
text('Parent',hS5E,'FontWeight','bold','FontSize',16,'String','o','units','centimeters','Position',[figSize.x - 3.8  figSize.y+0.4 0]);