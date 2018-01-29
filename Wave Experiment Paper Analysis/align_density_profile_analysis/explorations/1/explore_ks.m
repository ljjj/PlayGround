%%
% Simply use wave_speed obtained from previous analysis to shift profiles.
%% load data
dropboxFolder = 'C:\Users\jl2345\Dropbox (emonetlab)\';
manuscriptFolder = [dropboxFolder, 'users\setsu_kato\papers\xiongfei_setsu\manuscript\'];
load([manuscriptFolder, 'Figures\data for figs\data2.mat'])
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat'])
%%
plotProfiles = 1;
% delta_x = 122; % uncertainty in space due to image bin size when counting cells, um
% delta_t = 3.2; % uncertainty in time due to image taking frequency, s
Nexp = numel(data2);
zvec = -5000:100:8000;
rho = cell(1, Nexp);
for i = 20
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
    
    % obtained density pdf profiles
    rho{i} = zeros(N_s, numel(zvec)); 
    if plotProfiles
        figure('visible','off')
    end
    for j = 1:N_s % different sweeps
        % find moving coordinate
        z = cp(sweeps(j)).pos - wave_speed(i)*cp(sweeps(j)).time;
        
        % find back of wave
        cpc = cp(sweeps(j)).counts;
        ccpc = cumsum(cpc, 'reverse');
        back_ind = find(ccpc > total_cells_h(i), 1, 'last');
        
        
        % kernal smooth
        total = sum(cpc(back_ind:end));
        x = zeros(total, 1);
        x(1:ccpc(end)) = z(end);
        for k = numel(z):-1:(back_ind+1)
            x((ccpc(k)+1):ccpc(k-1)) = z(k-1);
        end
        [f, xi, bw] = ksdensity(x, zvec);
        rho{i}(j,:) = f;
        
        % show profiles
        if plotProfiles
            color = j/N_s * [1 0 0] + (N_s-j)/N_s * [0 0 1];
            nRow = 2;
            nCol = 3;
            
            subplot(nRow, nCol, 1)
            hold on
            plot(z, cpc/total_cells(i)/mean(diff(z)), 'color', color)
            xlabel('z')
            ylabel('cell counts')
            xlim([-5000 8000])
            ylim([0 1e-3])
            title('count')

            subplot(nRow, nCol, nCol+1)
            hold on
            plot(z, ccpc/total_cells(i), 'color', color)
            xlabel('z')
            ylabel('normalized cumulative cell counts')
            xlim([-5000 8000])
            ylim([0 1.2])

            subplot(nRow, nCol, 2)
            hold on
            plot(z, cp(sweeps(j)).counts_smooth/total_cells(i)/mean(diff(z)), 'color', color)
            xlabel('z')
            ylabel('cell counts')
            xlim([-5000 8000])
            ylim([0 1e-3])
            title('count smooth')

            subplot(nRow, nCol, nCol+2)
            hold on
            plot(z, cumsum(cp(sweeps(j)).counts_smooth, 'reverse')/total_cells(i), 'color', color)
            xlabel('z')
            ylabel('normalized cumulative cell counts')
            xlim([-5000 8000])
            ylim([0 1.2])

            subplot(nRow, nCol, 3)
            hold on
            plot(xi, f, 'color', color)
            xlabel('z')
            ylabel('pdf')
            xlim([-5000 8000])
            ylim([0 1e-3])
            title(['ksdensity bw = ', num2str(bw)])

            subplot(nRow, nCol, nCol+3)
            hold on
            plot(xi, 1-cumtrapz(xi, f), 'color', color)
            xlabel('z')
            ylabel('cdf')
            xlim([-5000 8000])
            ylim([0 1.2])
            end
    end
    if plotProfiles
        saveas(gcf, ['explore_ks_',num2str(i),'.png']);
        close(gcf)
    end
end