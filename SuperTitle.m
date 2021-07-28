function h = SuperTitle(titleText, titleArgs)
% puts one title at the top of a subplot
% Inputs: 
%   titleText = text input to title
%   titleArgs = cell of name-value pairs to pass into title()

if ~exist('titleArgs','var')
    titleArgs = {};
end

g = gca; % store current axes

set(gcf,'NextPlot','add');
axes; % make new

h = title(titleText, titleArgs{:}); % put title
set(gca,'Visible','off'); % hide current
set(h,'Visible','on'); % show title

axes(g); % restore previous axes
end