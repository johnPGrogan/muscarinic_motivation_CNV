function PlotTopoAndTrace(tVals,pVals,xTimes,varargin)
% function PlotTopoAndTrace(tVals,pVals,xTimes,varargin)
% 
% Plot a row of topoplots under a time series of those values
% 
% Inputs:
%   tVals = [nChans nTimes] matrix of t values (or beta, or voltage etc) to
%       use for plotting
%   pVals = [nChans nTimes] matrix of pValues
%   xTimes = [1 nTimes] vector of time points
%   varargin: name-param pairs for:
%       'topoTimes': timepoints to draw a topoplot for (default is 9 evenly
%           spaced across epoch)
%       'chansToPlot': cell array of channels to do time series for
%           (default is {'FPz','Cz','POz'})
%       'alpha': significance threshold (.05)
%       'xLab': time series x axis label ('time from ready cue (ms')
%       'ylab': time series y axis label ('t-statistic for incentive')
%       'doLegend': 1=draw legend, 0 = don't. (default =1)
%       'cbarArgs': cell array of args to pass to colorbar. default is
%           empty, which means do not draw a colorbar
%       'mapLimits': [min max] colorbar limits (10th and 90th percentiles)
%       'cmap': colormap (othercolour('PRGn11'))
%       'horizSpacing': horizontal spacing argument for subaxis (.01) -
%           requires 'subaxis.m', otherwise subplot.m is used
%       'sepFigs': whether to draw separate figures rather than both on one
%           subplot (default=0)
%       subplotInds = [rows, columns, index] for subplot/subaxis (default 
%           is [2 9 1])
%       title: title for line plot
% 
% John Grogan, 2021.


addpath ../../Experiments/MyFuncs/subAxis/ % use subaxis

load('chanlocs.mat');
%% parse
params = {  'topoTimes', linspace(min(xTimes), max(xTimes), 9);
            'chansToPlot', {'FPz', 'Cz', 'POz'};
            'alpha', .05;
            'xLab', 'time from ready cue (ms)';
            'yLab', 't-statistic for incentive';
            'doLegend',1;
            'cbarArgs',{};
            'mapLimits', [-1 1] .* max(abs(quantile(tVals, [.1 .9])),[],'all');
            'cmap', othercolor('PRGn11');
            'horizSpacing', .01;
            'sepFigs', 0;
            'subplotInds', [2 9 1];
            'title', '';
         };
     
[topoTimes, chansToPlot, alpha, xLab, yLab, doLegend, cbarArgs, mapLimits, cmap, horizSpacing, sepFigs, subplotInds,Title] = ...
    parsepvpairs(params(:,1), params(:,2), varargin{:});


%% plot a row of topoplots
if sepFigs
    figure();
end

for i = 1:length(topoTimes)
    if exist('subaxis','file')
        subaxis(subplotInds(1),subplotInds(2),subplotInds(3)*subplotInds(2)+i,'SpacingHorizontal', horizSpacing);
    else
        subplot(subplotInds(1),subplotInds(2),subplotInds(3)*subplotInds(2)+i);
    end
    topoplot(tVals(:,find(xTimes>=topoTimes(i),1,'first')), chanLocs(1:61),...
        'maplimits',mapLimits, 'colormap', cmap,...
        'style', 'both', 'plotrad',.55, 'headrad', .5,'emarker', {'.','k',[],1},...
        'numcontour', 6,'electrodes', 'on', 'nosedir', '+X',...
        'emarker2', {find(pVals(:,find(xTimes>=topoTimes(i),1,'first'))<alpha),'o','y',2});
    
%     title([num2str(round(xTimes(find(xTimes>=topoTimes(i),1,'first')))) ' ms']);
end

if ~isempty(cbarArgs)
    cbar = colorbar(cbarArgs{:}); % small
end



%% plot a few channels t-values

iChan = find(ismember(chanNames, chansToPlot)); % repeat it

if sepFigs
    figure();
else
    subaxis(subplotInds(1),subplotInds(2),(subplotInds(3)-1)*subplotInds(2)+1 : subplotInds(3)*subplotInds(2));
end

h = plot(xTimes, tVals(iChan,:),'LineWidth',3);
hold on;

xlim(round(minMax(xTimes,2)))
ylim(mapLimits);


yl = ylim;
j = 1;
for i = length(iChan):-1:1
    p = pbar(pVals(iChan(i),:),'alpha',alpha, 'xVals',xTimes, 'yVal', yl(1)+range(yl)*((j-1)/20), 'plotargs', {'LineWidth',5,'Color',h(i).Color});
    if any(~isnan(p.YData))
        j = j + 1;
    end    
end

yline(0, '--k','LineWidth',2);
xline(0, '--k','LineWidth',2);


xticks(round(topoTimes))

xlabel(xLab);
ylabel(yLab);%chanNames{iChan}]);

if doLegend
    legend(h, chanNames(iChan), 'Location','Best');
end

title(Title);
box off;
end