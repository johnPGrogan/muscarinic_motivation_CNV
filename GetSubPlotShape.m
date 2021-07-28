function [r, c] = GetSubPlotShape(n)
% function [r, c] = GetSubPlotShape(n)
% use ceil&sqrt to get good shape for subplots depending on number
% will be rectangular, roughly square
% 

r = ceil(sqrt(n));
c = ceil(n / r);

if nargout == 1 % return both in one number
    r = [r c];
end
end