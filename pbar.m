function h = pbar(pVals, varargin)
% h = pbar(pVals, varargin)
% plot a bar where p values are less than alpha
% Inputs: pVals is a vector of pvalues e.g. from tests at each time point
%         name-value pairs of optional inputs for graphical controls:
%           xVals = x-axis values to plot to (same size as pVals)
%           yVal = 0. y value to plot line along
%           alpha = .05. alpha threshold, only values below are plotted
%           plotargs = cell array of name-par to pass to plot (e.g. Color
%           or LineWidth)
% 
% Returns handle to plot
% 

names = {'yVal', 'alpha','plotargs','xVals'};
defaults = {0, .05, {'Color','k','LineWidth',5}, 1:length(pVals)};
[yVal, alpha, plotargs, xVals] = parsepvpairs(names, defaults, varargin{:});

pVal2 = double(pVals < alpha); % get less than alpha
pVal2(pVal2==0) = NaN; % set others to NaN

h = plot(xVals, pVal2 .* yVal, '-', plotargs{:});
end