function h = SuperimposeKS(data, ksArgs, yMinMax, colour, alpha)
% function SuperimposeKS(data, ksArgs, yMinMax, colour, alpha)
% 
% Superimpose a KSDensity of a distribution (e.g. RT) onto a figure.
% 
% Inputs:
%   data = data to run ksdensity() on. Will convert to column first
%   ksArgs = cell array of arguments passed into ksdensity as ksArgs{:},
%       e.g. {[0 100], 'Bandwidth',10};. default = {}
%   yMixMax = [min max] values of the filled object to appear on yaxis. It
%       will normalise the height to fill this at the peak. default = ylim
%   colour = colour to fill() (e.g. 'r' or [1 0 0]). default is from gca
%       colour order
%   alpha = alpha value of fill() (0-1). default = .5
% 
% Outputs:
%   h = handle to fill

if ~exist('ksArgs','var')
    ksArgs = {};
elseif ~iscell(ksArgs)
    error('ksArgs must be a cell array of arguments');
end
if ~exist('yMinMax','var')
    yMinMax = ylim();
elseif length(yMinMax)~=2
    error('yMinMax should be a 2-item vector');
end
if ~exist('colour','var')
    colour = get(gca,'ColorOrder');
    colour = colour(1,:);
end
if ~exist('alpha','var')
    alpha = .5;
end


[y,x] = ksdensity(data(:), ksArgs{:});

x1 = [x, fliplr(x)]; % mirror it
y1 = [y, zeros(size(y))]; % put zeros

% normalise y1 height
maxHeight = yMinMax(2) - yMinMax(1); % y units you want peak to reach
y1 = (y1 ./ max(y1,[],'all')) * maxHeight; % e.g. 0-2

% shift bottom up
y1 = y1 + yMinMax(1);

% plot
hold on;
h = fill(x1, y1, colour(1,:),'FaceAlpha',alpha);

end