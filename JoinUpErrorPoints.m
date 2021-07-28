function h1 = JoinUpErrorPoints(h, toJoin, plotargs)
% function h = JoinUpErrorPoints(h, toJoin)
% Join up pairs of errorbarplot points with lines
% 
% Inputs:
%   h: errorBarPlot handles for lines
%   toJoin: matrix of pairs of points to join up [n x 2]
%   plotargs: cell array of name-value pairs to pass into plot()
% 
% Outputs:
%   h1 = handles to new lines
% 

if ~exist('plotargs','var') || isempty(plotargs)
    plotargs = {};
end

n = size(toJoin,1);
hold on;
for i = 1:n
    for j = 1:numel(h)
        h1(j,i) = plot(h(j).XData(toJoin(i,:)), h(j).YData(toJoin(i,:)), ...
            'LineWidth', h(j).LineWidth, 'Color', h(j).Color, plotargs{:});
    end
end
