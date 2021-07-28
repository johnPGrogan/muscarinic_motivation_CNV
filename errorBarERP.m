function h = errorBarERP(EEG, elecToPlot, varargin)
% function h = errorBarERP(EEG, elecToPlot, varargin)
% use errorBarPlot on a subset of electrodes from EEG data (after epoching)
% Inputs:   EEG = eeglab structure after epoching
%           elecToPlot = cell array of channel names to plot
%           varargin: 'NAME', 'VALUE' pairs of optional arguments:
%                       'xLims' = xaxis limits, defaults to [min,max] of
%                           EEG.times
%                       'flipY' = Bool of whether to invert Y axis. detault
%                           is 0.
%                       'legend' = Bool of whether to make legend. Default
%                           is 0.
%                       'plotargs' = cell array of plotting arguments to be
%                           passed into errorBarPlot. Default is {}
% 
% Outputs:  h = errorBarPlot handle
% 
% Requires: errorBarPlot.m, parsepvpairs.m


% set deault param values unless optional values given
parNames = {'xLims','flipY','legend','plotargs'};
parVals  = {[min(EEG.times), max(EEG.times)], 0, 0, {}};
[xLims, flipY, isLegend, plotArgs] = parsepvpairs(parNames, parVals, varargin{:});

% get channels to plot
chanNames = {EEG.chanlocs(:).labels}; % names of electrodes
elec = ismember(chanNames,elecToPlot); % indices of electrodes to plot
if sum(elec) ~= length(elecToPlot)
    error('not all electrodes found in chanlocs')
end

% errorBarPlot those channels
plotData = permute(EEG.data(elec,:,:),[3,2,1]);
h = errorBarPlot(plotData,'area',1,'xaxisvalues',EEG.times,'plotargs',plotArgs);
xlim(xLims)

% plot zero x-axis
hold on
yline(0);

%plot zero y-axis (timelocking event)
yLims = ylim;
xline(0);

if flipY
    set(gca, 'YDir','reverse') % flip to plot negative up
end

if isLegend
    legend(h(:,1),{EEG.chanlocs(elec).labels},'Location','SouthWest')
end

end