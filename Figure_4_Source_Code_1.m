function PlotFigure4(allBetas, allPVals, xTimes2)

    if nargin==0
        % this one plots Figure 4
        load('Figure 4 - Source Data 1.mat', 'allBetas','allPVals','xTimes2')

        % these are for the figure supplements - see
        % TrioEEG_PlotAllFigures.m in the github for details
        % Figure 4 - Figure Supplement 1 - Source Data 1.mat
        % Figure 4 - Figure Supplement 2 - Source Data 1.mat
    end

    nVars = 3;
    
    load('chanlocs.mat');

    
    % plot
    topoTimes = [-100 300 650 1000 1350];
    timelockName = 'preparation cue'; % time from...
    
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
    for iV = 1:nVars
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