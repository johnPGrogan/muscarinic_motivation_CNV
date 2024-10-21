function PlotFigure3_S4(GNDEffects)
% Plot figure 3 - supplement 4
% 

    if nargin==0
        load('Figure 3 - Figure Supplement 4 - Source Data 1.mat','GNDEffects'); % cell array of GND outputs
    end


    topoTimes = [-100 300 650 1000 1350];
    timelockName = 'preparation cue onset'; % time from...


    % incentive
    if ~exist('othercolor','file')
        addpath '../../../../Documents/TCD/GeneralScripts/MatlabPackages/othercolor/';
    end
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
    