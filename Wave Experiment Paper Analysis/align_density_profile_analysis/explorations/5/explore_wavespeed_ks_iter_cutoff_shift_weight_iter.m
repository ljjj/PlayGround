%%
% Add a inner iteration to modify cutoff and shift mean until convergence,
% based on how close the individual KDE is to the mean, before regressing
% the wave speed with the outer iteration. Mixture model is removed.
% Convergence conditions are added for dynamic iterations.
%% load data
dropboxFolder = 'C:\Users\jl2345\Dropbox (emonetlab)\';
manuscriptFolder = [dropboxFolder, 'users\setsu_kato\papers\xiongfei_setsu\manuscript\'];
load([manuscriptFolder, 'Figures\data for figs\data2.mat'])
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat'])
%%
plotProfiles = 1;
bw = 200; % bandwidth of kernel smooth
nExp = numel(data2); % number of experiments
nMaxSpeedIter = 20; % number of iterations of wave speed
nMaxCutoffIter = 20; % number of iterations of cutoff
nShow = 1; % how many iterations to save figure
zVec = -5000:50:5000; % standard z vector to evaluate kernel smooth
nZ = numel(zVec);
for s = 2.^2 % hyperparameter, determines how strongly to downweight outlier pdf
dirName = [num2str(s)];
    if ~exist(dirName, 'dir')
        mkdir(dirName)
    end
    speeds = nan(nMaxSpeedIter+2, nExp); % record wave speed iterations
    speeds(1,:) = wave_speed;
    finalSpeed = nan(1, nExp); % final speed used to shift
    finalShift = nan(1, nExp); % constant shift within each experiment
    finalRho = nan(nZ, nExp); % weighted average of density pdf
    finalLabels = cell(1, nExp); % 1 = in the wave
    finalWeights = cell(1, nExp); % weight by closeness to mean
    finalCounts = cell(1, nExp); % cell counts
    finalWaveCounts = cell(1, nExp); % nubmer of cells in the wave
    
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
        sweepTimes = cellfun(@mean, {cp(sweeps).time});
        sampleCounts = vertcat(cp(sweeps).counts); % cell counts in each bin
        sampleX = vertcat(cp(sweeps).pos); % position of each bin
        
        % initialize iteration variables
        sweepWeights = ones(1, nSweeps)/nSweeps; % weight of each sweep
        sweepShifts = zeros(1, nSweeps); % constant shift after time shift c*t
        % find back of wave as label initialization
        initialSampleLabels = nan(nSweeps, nBins); % 1 = in the wave
        for j = 1:nSweeps
            cpc = cp(sweeps(j)).counts;
            ccpc = cumsum(cpc, 'reverse');
            backInd = find(ccpc > total_cells(i), 1, 'last');
            initialSampleLabels(j,1:backInd-1) = 0;
            initialSampleLabels(j,backInd:end) = 1;
        end
        if i == 2 % experiment 2 has cells not supposed to be at the end
            initialSampleLabels(:, cp(1).pos > 8000) = 0;
        end
        sampleLabels = initialSampleLabels; % probability to be in the wave sample
        
        % iterate speed and constant bias, mean pdf, weights, and mixture model
        sConverge = 0;
        for sIter = 0:(nMaxSpeedIter+1)
            
            % shift to moving frame
            sampleZ = sampleX - speeds(sIter+1,i)*vertcat(cp(sweeps).time) - sweepShifts(:);
            
            % kill extra cells in front of the wave
            if sIter > 0
                sampleLabels(sampleZ > 3000) = 0;
            end
            
            % aggregate samples
            weightedCounts = sweepWeights'.*sampleCounts; % cell counts are weighted by sweep weight
            waveWeights = weightedCounts.*sampleLabels; % sweep weighted cell counts in the wave component
            % kernal smooth wave sample in moving frame
            rhoAggregate = mvksdensity(sampleZ(:),zVec(:),'Bandwidth',bw,'Weights',waveWeights(:));
            
            % check speed convergence
            if sIter > 2
                sConverge = (speeds(sIter+1,i) - speeds(sIter,i))/speeds(sIter,i) < 0.001;
            end
            if sConverge
                countdown = countdown-1;
            else
                countdown = 2;
            end
            
            % only need to obtain estimated rho on standard vector in the last loop
            if sIter > nMaxSpeedIter || countdown == 0
                break
            end
            
            for cIter = 0:nMaxCutoffIter
                plotProfiles = mod(cIter,nShow) == 0 || descend == 0; % at regular iterations or when no gradient descent happens
                if plotProfiles
                    figure('visible','off')
                end
                
                % calculate wave shifts
                waveCounts = sampleLabels.*sampleCounts; % cell counts in wave component, sweep weight doesn't matter for each sweep is independent
                newSweepShifts = (sum(waveCounts.*sampleZ, 2)./sum(waveCounts, 2))';

                % KDE for each sweep
                sampleLabelsDeltaF = sampleLabels;
                sampleLabelsDeltaB = sampleLabels;
                rhoSweeps = nan(nZ, nSweeps); % calculated pdf estimated from inidividual wave sweeps
                rhoSweepsDeltaF = nan(nZ, nSweeps); % calculated pdf estimated from inidividual wave sweeps including 1 less bin
                rhoSweepsDeltaB = nan(nZ, nSweeps); % calculated pdf estimated from inidividual wave sweeps including 1 more bin
                cutoffInd = nan(1,nSweeps); % index where wave cuts off
                for j = 1:nSweeps
                    rhoSweeps(:,j) = mvksdensity(sampleZ(j,:)',zVec(:),'Bandwidth',bw,'Weights',waveCounts(j,:)');
                    cutoffInd(j) = find(sampleLabels(j,:) == 1, 1, 'first');
                    sampleLabelsDeltaF(j, cutoffInd(j)) = 0;
                    rhoSweepsDeltaF(:,j) = mvksdensity(sampleZ(j,:)',zVec(:),'Bandwidth',bw,'Weights',sampleLabelsDeltaF(j,:)'.*sampleCounts(j,:)');
                    if cutoffInd(j) > 1
                        sampleLabelsDeltaB(j, cutoffInd(j)-1) = 1;
                        rhoSweepsDeltaB(:,j) = mvksdensity(sampleZ(j,:)',zVec(:),'Bandwidth',bw,'Weights',sampleLabelsDeltaB(j,:)'.*sampleCounts(j,:)');
                    else
                        rhoSweepsDeltaB(:,j) = rhoSweeps(:,j);
                    end
                    
                    if plotProfiles
                        color = j/nSweeps * [1 0 0] + (nSweeps-j)/nSweeps * [0 0 1];

                        subplot(2,6,1)
                        hold on
                        plot(sampleX(j,:), sampleCounts(j,:), 'color', color)
                        plot(sampleX(j,:), sampleLabels(j,:)*10000, '--', 'color', color)
                        xlabel('x (\mum)')
                        ylabel('cell counts')
                        xlim([0 10000])

                        subplot(2,6,7)
                        hold on
                        plot(sampleX(j,:), cumsum(sampleCounts(j,:),'reverse'), 'color', color)
                        plot(sampleX(j,:), sampleLabels(j,:)*100000, '--', 'color', color)
                        xlabel('x (\mum)')
                        ylabel('cumulative cell counts')
                        xlim([0 10000])

                        subplot(2,6,2)
                        hold on
                        plot(sampleZ(j,:), sampleCounts(j,:), 'color', color)
                        plot(sampleZ(j,:), sampleLabels(j,:)*10000, '--', 'color', color)
                        xlabel('z (\mum)')
                        ylabel('cell counts (all)')
                        xlim([min(zVec) max(zVec)])

                        subplot(2,6,8)
                        hold on
                        plot(sampleZ(j,:), cumsum(sampleCounts(j,:),'reverse'), 'color', color)
                        plot(sampleZ(j,:), sampleLabels(j,:)*100000, '--', 'color', color)
                        xlabel('z (\mum)')
                        ylabel('cumulative cell counts (all)')
                        xlim([min(zVec) max(zVec)])

                        subplot(2,6,3)
                        hold on
                        plot(sampleZ(j,:), waveCounts(j,:), 'color', color)
                        xlabel('z (\mum)')
                        ylabel('cell counts (wave)')
                        xlim([min(zVec) max(zVec)])

                        subplot(2,6,9)
                        hold on
                        plot(sampleZ(j,:), cumsum(waveCounts(j,:),'reverse'), 'color', color)
                        xlabel('z (\mum)')
                        ylabel('cumulative cell counts (all)')
                        xlim([min(zVec) max(zVec)])

                        subplot(2,6,4)
                        hold on
                        plot(zVec, rhoSweeps(:,j), 'color', color)
                        xlabel('z (\mum)')
                        ylabel('ks pdf')
                        xlim([min(zVec) max(zVec)])
                        ylim([0 1e-3])

                        subplot(2,6,10)
                        hold on
                        plot(zVec, 1-cumtrapz(zVec, rhoSweeps(:,j)), 'color', color)
                        xlabel('z (\mum)')
                        ylabel('ks cdf')
                        xlim([min(zVec) max(zVec)])
                        ylim([0 1])
                    end
                end
                
                % calculate cutoff gradient
                L1Norms = sum(abs(rhoSweeps - rhoAggregate));
                L1NormsDeltaF = sum(abs(rhoSweepsDeltaF - rhoAggregate));
                L1NormsDeltaB = sum(abs(rhoSweepsDeltaB - rhoAggregate));
                descend = 0;
                for j = 1:nSweeps
                    if L1NormsDeltaF(j) < L1Norms(j)
                        sampleLabels(j, cutoffInd(j)) = 0;
                        descend = 1;
                    elseif L1NormsDeltaB(j) < L1Norms(j) && cutoffInd(j) > 1
                        sampleLabels(j, cutoffInd(j)-1) = 1;
                        descend = 1;
                    end
                end
                
                % calculate sweep weight
                epsilon = mean(L1Norms)/s;
                sweepWeights = exp(-L1Norms/epsilon);
                sweepWeights = sweepWeights / sum(sweepWeights);

                if plotProfiles
                    subplot(2,6,4)
                    hold on
                    plot(zVec, rhoAggregate, 'k--', 'lineWidth', 2)

                    subplot(2,6,10)
                    hold on
                    plot(zVec, 1-cumtrapz(zVec, rhoAggregate), 'k--', 'lineWidth', 2)

                    subplot(2,6,5)
                    plot(sweepTimes,sweepWeights)
                    xlabel('time (s)')
                    ylabel('weight')
                    xlim([0 ceil(max(sweepTimes)/100)*100])
                    ylim([0 0.2])
                    
                    subplot(2,6,6)
                    plot(sweepTimes,sum(waveCounts,2))
                    xlabel('time (s)')
                    ylabel('cell counts in wave')
                    xlim([0 ceil(max(sweepTimes)/100)*100])
                    yL = get(gca,'YLim');
                    yL(1) = 0;
                    ylim(yL)
                    
                    subplot(2,6,11)
                    plot(sweepTimes,sweepShifts + newSweepShifts)
                    xlabel('time (s)')
                    ylabel('shift (\mum)')
                    xlim([0 ceil(max(sweepTimes)/100)*100])
                    ylim([-1000 3000])

                    set(gcf, 'PaperPosition', [0 0 20 7]);
                    print(gcf, [dirName,'\exp_',num2str(i),'_sIter_',num2str(sIter),'_cIter_',num2str(cIter),'.png'], '-dpng', '-r300', '-opengl')
                    close(gcf)
                end
                
                if plotProfiles && ~ descend % finished plotting and descent
                    break
                end
            end
            
            % linear regress shifts on time to modify speed estimation
            sweepShifts = sweepShifts + newSweepShifts;
            x = sweepTimes;
            y = sweepShifts;
            meanX = sum(sweepWeights.*x);
            meanY = sum(sweepWeights.*y);
            meanXX = sum(sweepWeights.*x.^2);
            meanXY = sum(sweepWeights.*x.*y);
            beta = (meanXY - meanX*meanY)/(meanXX - meanX^2);
            sweepShifts(:) = meanY - beta * meanX;
            if sIter < nMaxSpeedIter+1
                speeds(sIter+2,i) = speeds(sIter+1,i) + beta;
            end
        end

        % finalize
        finalSpeed(i) = speeds(sIter+1,i);
        finalShift(i) = sweepShifts(1);
        finalRho(:,i) = rhoAggregate;
        finalLabels{i} = sampleLabels;
        finalWeights{i} = sweepWeights;
        finalCounts{i} = sampleCounts;
        finalWaveCounts{i} = sum(waveCounts,2);
    end
    save([dirName,'\results.mat'], 'bw', 'zVec', 'speeds', 'finalSpeed', 'finalShift', 'finalRho', 'finalLabels', 'finalWeights', 'finalCounts', 'finalWaveCounts');
end