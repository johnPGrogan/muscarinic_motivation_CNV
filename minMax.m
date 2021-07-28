function [Y, I] = minMax(X, dim)
% [Y, I] = minMax(X, dim)
% return minimum and maximum values
% Inputs:   X = vector or matrix of values
%           DIM = dimension to operate along for min() and max(), can be
%               integer or 'all'. see 'help min'
% 
% Outputs:  Y = multidimensional array of min and max values, concatenated
%               along DIM dimension, i.e. dim=1: 2xN array; dim=2: Nx2
%               array
%           I = indices returned from min() and max(). will not work if  
%               DIM = 'all'
% 
% John Grogan

if ~exist('dim','var') || isempty(dim) % default = 1
    dim = 1;
end

if isstr(dim) && strcmp(dim,'all') % cannot get inds for 'all'
    yMin = min(X, [], dim);    
    yMax = max(X, [], dim);
    I = [];
    
    % cat
    Y = cat(1, yMin, yMax);
else
    [yMin, iMin] = min(X, [], dim);
    [yMax, iMax] = max(X, [], dim);
    I = cat(dim, iMin, iMax);
    
    % cat
    Y = cat(dim, yMin, yMax);
end




