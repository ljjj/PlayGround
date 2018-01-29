%%
% Add a shift calculated from the mean of each wave_speed corrected pdf and
% regress with time to correct for wave speed. The regression weights are
% calculated based on L1 norm of each pdf with the weighted mean pdf.
%% load data
dropboxFolder = 'C:\Users\jl2345\Dropbox (emonetlab)\';
manuscriptFolder = [dropboxFolder, 'users\setsu_kato\papers\xiongfei_setsu\manuscript\'];
load([manuscriptFolder, 'Figures\data for figs\data2.mat'])
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat'])
%%
plotProfiles = 0;
% delta_x = 122; % uncertainty in space due to image bin size when counting cells, um
% delta_t = 3.2; % uncertainty in time due to image taking frequency, s
for scale = 2.^(1) % hyperparameter, determines how strongly to downweight outliers
if ~exist(['scale=',num2str(scale)], 'dir')
    mkdir(['scale=',num2str(scale)])
end
nExp = numel(data2); % number of experiments
nIter = 51; % number of iterations
nShow = 10; % how many iterations to save figure
zVec = -10000:100:10000;
finalSpeed = nan(1, nExp);
finalShift = nan(1, nExp);
finalRho = nan(nExp, numel(zVec));
for i = 1:nExp
    % initialization
    if isempty(data2(i).CellDensity)
        continue
    end
    cp = data2(i).CellDensity.CellPos;
    if isempty(cp(2).pos)
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    nS = numel(sweeps);
    
    % iterate to find the mean pdf as well as individual weights and shifts
    rhos = zeros(nS, numel(zVec));
    speed = wave_speed(i); % initialize wave speed estimates
    weights = ones(1, nS)/nS;
    shifts = zeros(1, nS);
    speeds = zeros(1, nIter);
    for iter = 1:nIter
        plotProfiles = mod(iter,nShow)==1;
        if plotProfiles
            figure('visible','off')
        end
        speeds(iter) = speed;
        % obtained density pdf profiles
        for j = 1:nS % different sweeps
            % find moving coordinate
            z = cp(sweeps(j)).pos - speed*cp(sweeps(j)).time - shifts(j);

            % find back of wave
            cpc = cp(sweeps(j)).counts;
            ccpc = cumsum(cpc, 'reverse');
            backInd = find(ccpc > total_cells_h(i), 1, 'last');

            % kernal smooth
            total = sum(cpc(backInd:end));
            x = zeros(total, 1);
            x(1:ccpc(end)) = z(end);
            for k = numel(z):-1:(backInd+1)
                x((ccpc(k)+1):ccpc(k-1)) = z(k-1);
            end
            [rhos(j,:), ~, bw] = ksdensity(x, zVec);
            
            % find mean shift
            shifts(j) = shifts(j) + trapz(zVec, rhos(j,:).*zVec);
            
            if plotProfiles
                color = j/nS * [1 0 0] + (nS-j)/nS * [0 0 1];
                
                subplot(2,2,1)
                hold on
                plot(zVec, rhos(j,:), 'color', color)
                xlabel('z (\mum)')
                ylabel('pdf')
                xlim([min(zVec) max(zVec)])
                ylim([0 1e-3])
                title(['ksdensity bw = ', num2str(bw)])

                subplot(2,2,3)
                hold on
                plot(zVec, 1-cumtrapz(zVec, rhos(j,:)), 'color', color)
                xlabel('z (\mum)')
                ylabel('cdf')
                xlim([min(zVec) max(zVec)])
                ylim([0 1.2])
            end
        end
        
        % calculate means and weights for next iteration
        meanRho = weights * rhos;
        L1Norms = sum(abs(rhos - meanRho),2)';
        epsilon = mean(L1Norms)/scale;
        weights = exp(-L1Norms/epsilon);
        weights = weights / sum(weights);
        x = cellfun(@mean, {cp(sweeps).time});
        y = shifts;
        
        if plotProfiles
            subplot(2,2,1)
            hold on
            plot(zVec, meanRho, 'k--')
            subplot(2,2,2)
            plot(x,weights)
            xlabel('time (s)')
            subplot(2,2,4)
            plot(x,y)        
            saveas(gcf, ['scale=',num2str(scale),'/explore_ks_wavespeed_iterate_scale=',num2str(scale),'_exp_',num2str(i),'_iter_',num2str(iter),'.png']);
            close(gcf)
        end
        
        % linear regress shifts on time to modify speed estimation
        meanX = sum(weights.*x);
        meanY = sum(weights.*y);
        meanXX = sum(weights.*x.^2);
        meanXY = sum(weights.*x.*y);
        beta = (meanXY - meanX*meanY)/(meanXX - meanX^2);
        speed = speed + beta;
        shifts(:) = meanY - beta * meanX;
    end
    
    % finalize
    finalSpeed(i) = speeds(end);
    finalShift(i) = shifts(1);
    for j = 1:nS % different sweeps
        % find moving coordinate
        z = cp(sweeps(j)).pos - finalSpeed(i)*cp(sweeps(j)).time - finalShift(i);

        % find back of wave
        cpc = cp(sweeps(j)).counts;
        ccpc = cumsum(cpc, 'reverse');
        backInd = find(ccpc > total_cells_h(i), 1, 'last');

        % kernal smooth
        total = sum(cpc(backInd:end));
        x = zeros(total, 1);
        x(1:ccpc(end)) = z(end);
        for k = numel(z):-1:(backInd+1)
            x((ccpc(k)+1):ccpc(k-1)) = z(k-1);
        end
        [rhos(j,:), ~, bw] = ksdensity(x, zVec);
    end
    % calculate means and weights for next iteration
    meanRho(i,:) = weights * rhos;
end
save(['scale=',num2str(scale),'\results.mat'], 'finalSpeed', 'finalShift', 'finalRho')
end