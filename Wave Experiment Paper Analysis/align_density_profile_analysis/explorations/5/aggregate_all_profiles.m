%% load data
dropboxFolder = 'C:\Users\jl2345\Dropbox (emonetlab)\';
manuscriptFolder = [dropboxFolder, 'users\setsu_kato\papers\xiongfei_setsu\manuscript\'];
load([manuscriptFolder, 'Figures\data for figs\data2.mat'])
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\wave_speed_analysis\final\wave_speeds_and_cell_counts.mat'])
s = 4;
load([manuscriptFolder, 'Figures\final figure May2016\fig2\fig2_material_May2017\align_density_profile_analysis\explorations\5\',num2str(s),'\results.mat'])
%% get moving frame data
nExp = numel(data2); % number of experiments
nZ = numel(zVec);

% find moving frame data
finalSampleZ = cell(1, nExp);
for i = 1:nExp
    % initialize
    if isempty(data2(i).CellDensity)
        continue
    end
    cp = data2(i).CellDensity.CellPos;
    if isempty(cp(2).pos) % for some experiments data is in every other struct
        sweeps = 1:2:numel(cp);
    else
        sweeps = 1:numel(cp);
    end
    
    % shift to moving frame
    finalSampleZ{i} = vertcat(cp(sweeps).pos) - finalSpeed(i)*vertcat(cp(sweeps).time) - finalShift(i);
end
%% aggregate data for each condition and iterate shift
nConditions = 3;
rhoAggregateShifted = nan(nConditions, nZ);
expSelect = {[1 3 14 15 16] [4 5 6 17 18] [7 8 9 19 20]};
dz = zVec(2)-zVec(1);
nMaxIter = 20;
deltaZ = zeros(1,nExp);

for c = 1:nConditions
    exps = expSelect{c};
    nE = numel(exps);
    
    % aggregate data
    sampleZ = vertcat(finalSampleZ{exps});
    sampleLabels = vertcat(finalLabels{exps});
    sweepWeights = horzcat(finalLabels{exps});
    sampleCounts = vertcat(finalCounts{exps});
    waveWeights = sweepWeights'.*sampleCounts.*sampleLabels; % sweep weighted cell counts in the wave component
    zi = 0;
    zInds = cell(1,nE); % index in aggregate data corresponding to each experiment
    for e = 1:nE
        nSweeps = size((finalSampleZ{exps(e)}),1);
        zInds{e} = zi+1:zi+nSweeps;
        zi = zi + nSweeps;
    end
    
    % iterate to adjust shift
    rhoExps = finalRho(:,exps); % KDE estimated from inidividual experiments
    for iter = 0:nMaxIter
        % find aggregate density
        rhoAggregateShifted(c,:) = mvksdensity(sampleZ(:),zVec(:),'Bandwidth',bw,'Weights',waveWeights(:));
    
        % get KDE for shifted individual exp
        rhoExpsDeltaF = [zeros(1,nE); rhoExps(1:end-1,:)]; % pdf estimated from inidividual experiments shifted forward
        rhoExpsDeltaB = [rhoExps(2:end,:); zeros(1,nE)]; % pdf estimated from inidividual experiments shifted backward

        % calculate shift gradient
        L1Norms = sum(abs(rhoExps - rhoAggregateShifted(c,:)'));
        L1NormsDeltaF = sum(abs(rhoExpsDeltaF - rhoAggregateShifted(c,:)'));
        L1NormsDeltaB = sum(abs(rhoExpsDeltaB - rhoAggregateShifted(c,:)'));
        descend = 0;
        for e = 1:nE
            if L1NormsDeltaF(e) < L1Norms(e)
                sampleZ(zInds{e},:) = sampleZ(zInds{e},:) + dz;
                deltaZ(exps(e)) = deltaZ(exps(e)) + dz;
                rhoExps(:,e) = rhoExpsDeltaF(:,e);
                descend = 1;
            elseif L1NormsDeltaB(e) < L1Norms(e)
                sampleZ(zInds{e},:) = sampleZ(zInds{e},:) - dz;
                deltaZ(exps(e)) = deltaZ(exps(e)) - dz;
                rhoExps(:,e) = rhoExpsDeltaB(:,e);
                descend = 1;
            end
        end
        
        clrs = {'r', [1 0.843137254901961 0], 'b'};

        figure('visible','off')
        subplot(2,1,1)
        hold on
        for e = 1:nE
            plot(zVec, rhoExps(:,e), 'color', clrs{c})
        end
        plot(zVec, rhoAggregateShifted(c,:), 'k--', 'lineWidth', 2)
        xlabel('z (\mum)')
        ylabel('pdf')
        xlim([-2500 2000])

        subplot(2,1,2)
        hold on
        for e = 1:nE
            plot(zVec, cumsum(rhoExps(:,e),'reverse')*dz, 'color', clrs{c})
        end
        plot(zVec, cumsum(rhoAggregateShifted(c,:),'reverse')*dz, 'k--', 'lineWidth', 2)
        xlabel('z (\mum)')
        ylabel('pdf')
        xlim([-2500 2000])

        set(gcf, 'PaperPosition', [0 0 10 6]);
        print(gcf, [num2str(s),'\aggregate\exp_density_profiles_cond_',num2str(c),'_iter_',num2str(iter),'.png'], '-dpng', '-r300', '-opengl')
        
        if ~descend
            break
        end
    end
end
%% aggregate data from experiments for final profile and bootstrapped epdf
binEdges = [(zVec-dz/2) zVec(end)+dz/2];
nBootstrap = 10000;
bsPDF = nan(3, nBootstrap, nZ);

for c = 1:nConditions
    exps = expSelect{c};
    
    % aggregate data
    sampleZ = vertcat(finalSampleZ{exps});
    sampleLabels = vertcat(finalLabels{exps});
    sweepWeights = horzcat(finalLabels{exps});
    sampleCounts = vertcat(finalCounts{exps});
    waveWeights = sweepWeights'.*sampleCounts.*sampleLabels; % sweep weighted cell counts in the wave component
    
    % find aggregate density
    rhoAggregateShifted(c,:) = mvksdensity(sampleZ(:),zVec(:),'Bandwidth',bw,'Weights',waveWeights(:));

    % bootstrap epdf
%     bsPDF(c,:,:) = bootstrp(nBootstrap, @(x) wepdf(x, binEdges), sampleZ(:), 'Weights', waveWeights(:));
    bsPDF(c,:,:) = bootstrp(nBootstrap, @(x) mvksdensity(x,zVec(:),'Bandwidth',bw), sampleZ(:), 'Weights', waveWeights(:));
end
bsCDF = cumsum(bsPDF, 3, 'reverse') * dz;
%% calculate CDF error from bootstrap
bsPDFMean = nan(nConditions, nZ);
bsPDFStd = nan(nConditions, nZ);

bsCDFMean = nan(nConditions, nZ);
bsCDFStd = nan(nConditions, nZ);

for c = 1:nConditions
    bsPDFMean(c,:) = mean(bsPDF(c,:,:));
    bsPDFStd(c,:) = std(bsPDF(c,:,:));
    bsCDFMean(c,:) = mean(bsCDF(c,:,:));
    bsCDFStd(c,:) = std(bsCDF(c,:,:));
end

bsCDFStdMax = nan(1,nConditions);
for c = 1:nConditions
    bsCDFStdMax(c) = max(bsCDFStd(c,:));
    disp(['Max std in CDF of condition ', num2str(c), ' is ', num2str(bsCDFStdMax(c))])
end

bsCDFAbsDiffMax = nan(nConditions, nConditions);
for c1 = 1:nConditions
    for c2 = 1:nConditions
        bsCDFAbsDiffMax(c1, c2) = max(abs(bsCDFMean(c1,:) - bsCDFMean(c2,:)));
        disp(['Max abs diff in CDF of condition ', num2str(c1), ' and condition ',...
            num2str(c2), ' is ', num2str(bsCDFAbsDiffMax(c1, c2))])
    end
end
%% save
save([num2str(s),'\aggregate\bootstrappedData.mat'], 'rhoAggregateShifted', 'deltaZ', 'nConditions', 'zVec', 'bsPDF', 'bsCDF', 'bsPDFMean', 'bsPDFStd', 'bsCDFMean', 'bsCDFStd', 'bsCDFStdMax', 'bsCDFAbsDiffMax')
%% load
load([num2str(s),'\aggregate\bootstrappedData.mat'])
%% show PDFs and CDFs
clrs = {'r', [1 0.843137254901961 0], 'b'};

figure
subplot(2,1,1)
hold on
for c = 1:nConditions
    plot(zVec, rhoAggregateShifted(c,:), 'color', clrs{c})
    plot(zVec, bsPDFMean(c,:), '--', 'color', clrs{c})
    patch([zVec flip(zVec)],[bsPDFMean(c,:)+bsPDFStd(c,:) flip(bsPDFMean(c,:)-bsPDFStd(c,:))], clrs{c},'EdgeColor','none','LineWidth',1,'EdgeAlpha',0.3,'FaceAlpha',0.3);
end
xlabel('z (\mum)')
ylabel('pdf')
xlim([-2500 2000])

subplot(2,1,2)
hold on
for c = 1:nConditions
    plot(zVec, cumsum(rhoAggregateShifted(c,:),'reverse')*dz, 'color', clrs{c})
    plot(zVec, bsCDFMean(c,:), '--', 'color', clrs{c})
    patch([zVec flip(zVec)],[bsCDFMean(c,:)+bsCDFStd(c,:) flip(bsCDFMean(c,:)-bsCDFStd(c,:))], clrs{c},'EdgeColor','none','LineWidth',1,'EdgeAlpha',0.3,'FaceAlpha',0.3);
end
xlabel('z (\mum)')
ylabel('pdf')
xlim([-2500 2000])

set(gcf, 'PaperPosition', [0 0 10 6]);
print(gcf, [num2str(s),'\aggregate\aggregate_density_profiles.png'], '-dpng', '-r300', '-opengl')
%% show CDF differences colorbar
figure
b = bar3(bsCDFAbsDiffMax);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
xlabel('Condition #')
ylabel('Condition #')
caxis([0 0.15])
zlim([0 0.15])
title('sup|CDF(i)-CDF(j)|')

set(gcf, 'PaperPosition', [0 0 6 6]);
print(gcf, [num2str(s),'\aggregate\aggregate_KS.png'], '-dpng', '-r300', '-opengl')
%% calculate CDF error from kernel smooth of each experiment
ksCDF = cumsum(finalRho, 'reverse') * dz;

ksCDFMean = nan(nExp, nZ);
ksCDFStd = nan(nExp, nZ);

for e = 1:nExp
    ksCDFMean(e,:) = nanmean(ksCDF(e,:,:));
    ksCDFStd(e,:) = nanstd(ksCDF(e,:,:));
end

ksCDFStdMax = nan(1,nExp);
for e = 1:nExp
    ksCDFStdMax(e) = max(ksCDFStd(e,:));
    disp(['Max std in CDF of condition ', num2str(e), ' is ', num2str(ksCDFStdMax(e))])
end

ksCDFAbsDiffMax = nan(nExp, nExp);
for c1 = 1:nExp
    for c2 = 1:nExp
        ksCDFAbsDiffMax(c1, c2) = max(abs(ksCDFMean(c1,:) - ksCDFMean(c2,:)));
        disp(['Max abs diff in CDF of condition ', num2str(c1), ' and condition ',...
            num2str(c2), ' is ', num2str(ksCDFAbsDiffMax(c1, c2))])
    end
end
%% show profiles
clrs = {'r', [1 0.843137254901961 0], 'b'};

figure
subplot(2,1,1)
hold on
for c = 1:nConditions
    for e = expSelect{c}
        plot(zVec, finalRho(:,e), 'color', clrs{c})
    end
end
xlabel('z (\mum)')
ylabel('pdf')
xlim([-2500 2000])

subplot(2,1,2)
hold on
for c = 1:nConditions
    for e = expSelect{c}
        plot(zVec, ksCDF(:,e), 'color', clrs{c})
    end
end
xlabel('z (\mum)')
ylabel('pdf')
xlim([-2500 2000])

set(gcf, 'PaperPosition', [0 0 10 6]);
print(gcf, [num2str(s),'\aggregate\exp_density_profiles.png'], '-dpng', '-r300', '-opengl')
%% show CDF differences colorbar
figure
b = bar3(ksCDFAbsDiffMax);
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
xlabel('Exp #')
ylabel('Exp #')
caxis([0 1e-8])
zlim([0 1e-8])
title('sup|CDF(i)-CDF(j)|')

set(gcf, 'PaperPosition', [0 0 6 6]);
print(gcf, [num2str(s),'\aggregate\exp_KS.png'], '-dpng', '-r300', '-opengl')