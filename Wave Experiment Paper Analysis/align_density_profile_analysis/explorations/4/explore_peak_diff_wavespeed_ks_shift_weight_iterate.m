%%
% Tried to find different wave speed for each sweep but peak identificaiton
% is too unreliable to do this.
%% load data
dropboxFolder = 'C:\Users\jl2345\Dropbox (emonetlab)\';
manuscriptFolder = [dropboxFolder, 'users\setsu_kato\papers\xiongfei_setsu\manuscript\'];
load([manuscriptFolder, 'Figures\data for figs\data2.mat'])
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat'])
%%
plotProfiles = 1;
bw = 200; % bandwidth of kernel smooth
nExp = numel(data2); % number of experiments
nIter = 10; % number of iterations
nShow = 1; % how many iterations to save figure
zVec = -10000:100:10000; % standard z vector to evaluate kernel smooth
nZ = numel(zVec);
for s = 2.^1 % hyperparameter, determines how strongly to downweight outlier pdf
dirName = num2str(s);
    if ~exist(dirName, 'dir')
        mkdir(dirName)
    end
    finalSpeed = nan(1, nExp); % final speed used to shift
    finalShift = nan(1, nExp); % constant shift within each experiment
    finalRho = nan(nZ, nExp); % weighted average of density pdf
    
    for i = 1:nExp
        % initialize constants
        if isempty(data2(i).CellDensity)
            continue
        end
        cp = data2(i).CellDensity.CellPos;
        if isempty(cp(2).pos) % for some experiments data is in every other struct
            sweeps = 1:2:numel(cp);
        else
            sweeps = 1:numel(cp);
        end
        nSweeps = numel(sweeps); % number of sweeps
        nBins = numel(cp(1).pos); % number of cell count bins in each sweep
        sampleCounts = nan(nSweeps-2, nBins); % cell counts in each bin
        for j = 1:nSweeps-2
            sampleCounts(j,:) = cp(sweeps(j+1)).counts;
        end
        sampleX = vertcat(cp(sweeps(2:end-1)).pos); % position of each bin
        
        % initialize iteration variables
        sweepWeights = ones(1, nSweeps-2)/(nSweeps-2); % weight of each sweep
        sweepShifts = zeros(1, nSweeps-2); % constant shift after time shift c*t
        
        % find peak wave speed
        peakPos = zeros(1,nSweeps);
        peakTime = zeros(1,nSweeps);
        for j = 1:nSweeps
            cm = cp(j).counts_smooth; % later wave for accuracy but not too late to keep more cells
            dsd = diff(sign( diff(cm) ));
            den_max_ind = find(dsd == -2) + 1; % find where first derivative = 0 and second derivative < 0
            if isempty(den_max_ind)
                den_max_ind = 1;
            end
            den_max_ind = den_max_ind(end);
            peakPos(j) = cp(sweeps(j)).pos(den_max_ind);
            peakTime(j) = cp(sweeps(j)).time(den_max_ind);
        end
        peakSpeed = (peakPos(3:end) - peakPos(1:end-2))./(peakTime(3:end) - peakTime(1:end-2));
        
        % iterate speed and constant bias, mean pdf, weights, and mixture model
        for iter = 0:(nIter+1)
            plotProfiles = mod(iter,nShow)==0;
            if plotProfiles
                figure('visible','off')
            end
            
            % shift to moving frame
            sampleZ = sampleX - peakSpeed(:)*vertcat(cp(sweeps(2:end-1)).time) - sweepShifts(:);
            
            % find back of wave as new label
            sampleLabels = nan(nSweeps-2, nBins); % 1 = in the wave
            for j = 1:nSweeps-2
                cpc = cp(sweeps(j+1)).counts;
                ccpc = cumsum(cpc, 'reverse');
                backInd = find(ccpc > total_cells_h(i), 1, 'last');
                initialSampleLabels(j,1:backInd-1) = 0;
                initialSampleLabels(j,backInd:end) = 1;
            end
            
            weightedCounts = sweepWeights'.*sampleCounts; % cell counts are weighted by sweep weight
            waveWeights = weightedCounts.*sampleLabels; % sweep weighted cell counts in the wave component
            % kernal smooth wave sample in moving frame
            rhoAggregate = mvksdensity(sampleZ(:),zVec(:),'Bandwidth',bw,'Weights',waveWeights(:));
            if iter > nIter % only need to obtain estimated rho on standard vector in the last loop
                break
            end
            
            
            % calculate wave shifts based on new label
            waveCounts = sampleLabels.*sampleCounts; % cell counts in wave component, sweep weight doesn't matter for each sweep is independent
            sweepShifts(:) = sweepShifts(:) + sum(waveCounts.*sampleZ, 2)./sum(waveCounts, 2);

            % calculate sweep weights for next iteration
            rhoSweeps = nan(nZ, nSweeps); % calcuated pdf estimated from inidividual sweeps
            for j = 1:nSweeps
                rhoSweeps(:,j) = mvksdensity(sampleZ(j,:)',zVec(:),'Bandwidth',bw,'Weights',sampleLabels(j,:)'.*sampleCounts(j,:)');
                
                if plotProfiles
                    color = j/nSweeps * [1 0 0] + (nSweeps-j)/nSweeps * [0 0 1];
                    
                    subplot(2,5,1)
                    hold on
                    plot(sampleX(j,:), sampleCounts(j,:), 'color', color)
                    plot(sampleX(j,:), sampleLabels(j,:)*10000, '--', 'color', color)
                    xlabel('x (\mum)')
                    ylabel('cell counts')
                    xlim([0 10000])

                    subplot(2,5,6)
                    hold on
                    plot(sampleX(j,:), cumsum(sampleCounts(j,:),'reverse'), 'color', color)
                    plot(sampleX(j,:), sampleLabels(j,:)*100000, '--', 'color', color)
                    xlabel('x (\mum)')
                    ylabel('cumulative cell counts')
                    xlim([0 10000])
                    
                    subplot(2,5,2)
                    hold on
                    plot(sampleZ(j,:), sampleCounts(j,:), 'color', color)
                    plot(sampleZ(j,:), sampleLabels(j,:)*10000, '--', 'color', color)
                    xlabel('z (\mum)')
                    ylabel('cell counts (all)')
                    xlim([min(zVec) max(zVec)])

                    subplot(2,5,7)
                    hold on
                    plot(sampleZ(j,:), cumsum(sampleCounts(j,:),'reverse'), 'color', color)
                    plot(sampleZ(j,:), sampleLabels(j,:)*100000, '--', 'color', color)
                    xlabel('z (\mum)')
                    ylabel('cumulative cell counts (all)')
                    xlim([min(zVec) max(zVec)])
                    
                    subplot(2,5,3)
                    hold on
                    plot(sampleZ(j,:), waveCounts(j,:), 'color', color)
                    xlabel('z (\mum)')
                    ylabel('cell counts (wave)')
                    xlim([min(zVec) max(zVec)])
                    
                    subplot(2,5,8)
                    hold on
                    plot(sampleZ(j,:), cumsum(waveCounts(j,:),'reverse'), 'color', color)
                    xlabel('z (\mum)')
                    ylabel('cumulative cell counts (all)')
                    xlim([min(zVec) max(zVec)])
                    
                    subplot(2,5,4)
                    hold on
                    plot(zVec, rhoSweeps(:,j), 'color', color)
                    xlabel('z (\mum)')
                    ylabel('ks pdf')
                    xlim([min(zVec) max(zVec)])
                    ylim([0 1e-3])

                    subplot(2,5,9)
                    hold on
                    plot(zVec, 1-cumtrapz(zVec, rhoSweeps(:,j)), 'color', color)
                    xlabel('z (\mum)')
                    ylabel('ks cdf')
                    xlim([min(zVec) max(zVec)])
                    ylim([0 1])
                end
            end
            L1Norms = sum(abs(rhoSweeps - rhoAggregate));
            epsilon = mean(L1Norms)/s;
            sweepWeights = exp(-L1Norms/epsilon);
            sweepWeights = sweepWeights / sum(sweepWeights);
            x = cellfun(@mean, {cp(sweeps).time});
            y = sweepShifts;
            % linear regress shifts on time to modify speed estimation
            meanX = sum(sweepWeights.*x);
            meanY = sum(sweepWeights.*y);
            meanXX = sum(sweepWeights.*x.^2);
            meanXY = sum(sweepWeights.*x.*y);
            beta = (meanXY - meanX*meanY)/(meanXX - meanX^2);
            sweepShifts(:) = meanY - beta * meanX;
            if iter < nIter+1
                speeds(iter+2,i) = speeds(iter+1,i) + beta;
            end
            
            if plotProfiles
                subplot(2,5,4)
                hold on
                plot(zVec, rhoAggregate, 'k--', 'lineWidth', 2)
                
                subplot(2,5,9)
                hold on
                plot(zVec, 1-cumtrapz(zVec, rhoAggregate), 'k--', 'lineWidth', 2)
                
                subplot(2,5,5)
                plot(x,sweepWeights)
                xlabel('time (s)')
                ylabel('weight')
                xlim([0 ceil(max(x)/100)*100])
                ylim([0 0.2])
                
                subplot(2,5,10)
                plot(x,y)
                xlabel('time (s)')
                ylabel('shift (\mum)')
                xlim([0 ceil(max(x)/100)*100])
                ylim([-1000 3000])
                
                set(gcf, 'PaperPosition', [0 0 20 7]);
                print(gcf, [dirName,'\exp_',num2str(i),'_iter_',num2str(iter),'.png'], '-dpng', '-r300', '-opengl')
                close(gcf)
            end
        end

        % finalize
        finalSpeed(i) = speeds(end,i);
        finalShift(i) = sweepShifts(1);
        finalRho(:,i) = rhoAggregate;
    end
    save([dirName,'\results.mat'], 'speeds', 'finalSpeed', 'finalShift', 'finalRho');
end