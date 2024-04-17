%% plot clusterperm fig5
% load up combined files with permutations
% cluster to get p-vals
% plot

clear; 
%% combine
if ~exist('spatial_neighbors','file'); eeglab nogui; end

timelock = 'readyBS';
file1 = sprintf('RegressPermuteStuff2_%s.mat',timelock);

if ~exist('xTimes','var'); load(file1, 'xTimes','timesToPlot'); end

saveFolder = './Data';
% saveFolder = 'C:\Users\groganj1\OneDrive - TCDUD.onmicrosoft.com\Oxford\AnalysisOD\Trio_EEG\Data';

vNames = {'velres', 'srt','depAngle'};
nVars = length(vNames);
% CombinePermData(saveFolder,vNames) % makes _combined files


%%

load('chanlocs.mat','chanLocs');

nT = 410;
nChans = 61;

for iV = 1:nVars
    if exist(sprintf('%s_regressCluster3_%s.mat',timelock, vNames{iV}),'file');
%         continue; 
    end

    f = what(saveFolder);
    files = f.mat(cellRegexpi(f.mat, sprintf('WholeBrainRegressPermuteParallel3_combined_%s',vNames{iV}))>0);
    
    % this should now just cat them, across nperm sets
    d.tValsPerm = [];
    for i = 1:length(files)
        data = load(fullfile(saveFolder, files{i}));
        d.tValsPerm = cat(3, d.tValsPerm, data.tValsPerm);
    end

    % now get real values

    filesReal = f.mat(cellRegexpi(f.mat, sprintf('WholeBrainRegressReal2_%s', vNames{iV}))>0);
    data = load(fullfile(saveFolder, filesReal{1}));
    d.tVals = data.tVals;
    d.betaVals = data.betaVals;

    if size(d.tVals,1) > nT
        t = xTimes(isBetween(xTimes, data.timesToPlot)); % get true-reg times
        d.tVals(~isBetween(t,timesToPlot),:,:) = []; % remove outliers
        d.betaVals(~isBetween(t,timesToPlot),:,:) = []; % remove outliers
    end
    
  

    %% or could I do clustering? using funcs from Mass_Univariate_Toolbox
    df = 18576; % degrees of freedom in glme - or should this be nPP-1?
    tail = 0;
    tic;
    d.pVals = FindClustersLikeGND(d.tVals, d.tValsPerm, chanLocs, tail, df);
    disp(toc/60);
    
    %% save this
    d = rmfield(d, 'tValsPerm');
    save(sprintf('%s_regressCluster3_%s.mat',timelock, vNames{iV}), '-struct','d');
end


%% now load

%% now load the regressions of each var against reward
% % timelock = 'targ';
% load(sprintf('trioEEGAnalyseByTrials_%s.mat',timelock),'xTimes', 'timesToPlot')
% % timesToPlot = [-200 1500];
xTimes2 = xTimes(isBetween(xTimes, timesToPlot));

% vNames = {'velres','srt','depAngle'};%,'endptVarCorr'};

[allBetas, allPVals, allTVals] = deal(zeros(sum(isBetween(xTimes,timesToPlot)),61,length(vNames)));
for iV = 1:nVars
    a = load(sprintf('%s_regressCluster3_%s.mat',timelock, vNames{iV}));
    allBetas(:,:,iV) = a.betaVals(1:410,:);
    allTVals(:,:,iV) = a.tVals(1:410,:);
    allPVals(:,:,iV) = a.pVals(1:410,:);    
end


%% plot
timelock = 'readyBS';
topoTimes = [-100 300 650 1000 1350];
timelockName = 'saccade onset'; % time from...

if ~exist('othercolor','file')
    addpath '../../../../Documents/TCD/GeneralScripts/MatlabPackages/othercolor/';
end
cmap = othercolor('PRGn11');

alpha = .05 ./ numel(allPVals(:,:,:,1));

mapLims = [-.15 .15] ;
chansToPlot = {'FPz','Cz','POz'};
yLabs = {'residual peak velocity','RT','distractor pull','endptVarCorr'};

% figure();
% close all;
set(0, 'DefaultAxesFontSize',14);
for iV = 1%:nVars
    if iV==1
        cbarArgs = {'Position',[0.0625 0.1952 0.0232 0.1905]};
    else
        cbarArgs = {};
    end
    f = figure();
    set(f, 'DefaultAxesFontSize',20);

    

    PlotTopoAndTrace(allBetas(:,:,iV)', allPVals(:,:,iV)',xTimes2,...
        'topoTimes', topoTimes, 'chansToPlot',chansToPlot,...
        'alpha', alpha, 'mapLimits',mapLims, 'cmap', cmap,...
        'yLab', '\beta coefficient', 'title', yLabs{iV},...
        'doLegend',iV==1, 'cbarArgs',cbarArgs, 'subplotInds',[2 length(topoTimes) 1],...
        'xLab', sprintf('time from %s (ms)', timelockName), 'sepFigs',0);
    box off;

    % save figures
    for j = 1:length(f)
%         saveas((f(j)), sprintf('./Figs/PlotClusterPerms3_%s_%s_%d.svg', timelock, varNames{iV}, j))
    end
end

