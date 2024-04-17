
clc; clear; eeglab; close all;

set(0, 'DefaultAxesFontSize',20)
%% do the cluster-permutation effects
timelock = 'readyBS';
topoTimes = [-100 300 650 1000 1350];
timelockName = 'preparation cue onset'; % time from...

load(sprintf('TrioEEGAnalyseByTrials_%s.mat',timelock), 'GNDEffects');

    %% incentive
    cmap = othercolor('PRGn11');
    mapLims = [-7.5 7.5];
    chansToPlot = {'FPz','Cz','POz'};
    yLabs = {'incentive', [], 'THP', [], 'THP * incentive'};
    

    for i = [1 3 5]
        if i==1
            cbarArgs = {'Position',[0.0625 0.1952 0.0232 0.1905]};
        else
            cbarArgs = {};
        end
        f = figure(); 
        set(f, 'DefaultAxesFontSize',14);
        PlotTopoAndTrace(GNDEffects{i}.grands_t(:,:,end), GNDEffects{i}.perm_t_tests(end).adj_pval,...
            GNDEffects{i}.time_pts,...
            'topoTimes', topoTimes, 'chansToPlot',chansToPlot, 'mapLimits',mapLims,...
            'yLab', 't-statistic', 'title', yLabs{i}, 'cmap', cmap,...
            'doLegend', i==1, 'cbarArgs',cbarArgs, 'subplotInds', [2 length(topoTimes) 1],...
            'sepFigs',0,'xLab', sprintf('time from %s (ms)', timelockName));
    end
    


%% now load the regressions of each var against reward
% timelock = 'targ';
load(sprintf('trioEEGAnalyseByTrials_%s.mat',timelock),'xTimes', 'timesToPlot')
% timesToPlot = [-200 1500];
xTimes2 = xTimes(isBetween(xTimes, timesToPlot));

vNames = {'velres','srt','depAngle'};%,'endptVarCorr'};

[allBetas, allPVals, allTVals] = deal(zeros(sum(isBetween(xTimes,timesToPlot)),61,length(vNames)));
for i = 1:length(vNames)
    a = load(sprintf('%s_regressCluster_%s.mat',timelock, vNames{i}));
    allBetas(:,:,i) = a.betaVals(1:410,:);
    allTVals(:,:,i) = a.tVals(1:410,:);
    allPVals(:,:,i) = a.pVals(1:410,:);    
end


%%
cmap = othercolor('PRGn11');

alpha = .05 ./ numel(allPVals(:,:,:,1));

mapLims = [-.15 .15] ;
chansToPlot = {'FPz','Cz','POz'};
yLabs = {'residual peak velocity','RT','distractor pull','endptVarCorr'};

figure();
for i = 1:3%length(vNames)
    if i==1
        cbarArgs = {'Position',[0.0625 0.1952 0.0232 0.1905]};
    else
        cbarArgs = {};
    end
    f = figure();
    set(f, 'DefaultAxesFontSize',14);

    PlotTopoAndTrace(allBetas(:,:,i)', allPVals(:,:,i)',xTimes2,...
        'topoTimes', topoTimes, 'chansToPlot',chansToPlot,...
        'alpha', alpha, 'mapLimits',mapLims, 'cmap', cmap,...
        'yLab', '\beta coefficient', 'title', yLabs{i},...
        'doLegend',i==1, 'cbarArgs',cbarArgs, 'subplotInds',[2 length(topoTimes) 1],...
        'xLab', sprintf('time from %s (ms)', timelockName));
    box off;
end


%% do an empty topoplot with legend colours

f = figure();
load('chanlocs.mat');
inds = find(ismember(chanNames,chansToPlot));
cols = get(gca,'ColorOrder');
subplot(2,5,6);
for i = 1:length(inds)
    topoplot(zeros(61,1), chanLocs(1:61),'mapLimits',[-1 1], 'colormap', crameri('vik'), 'emarker2', {inds(i), 'o', cols(i,:), 6});
    hold on;
end
