function yes = isBetween(data, limits, isInclusive)
% check whether number lies between two others
% Inputs:
%   data = matrix of numbers
%   limits = [min max] limits to look between
%   isInclusive: 1 = inclusive boundaries, 0 = exclusive. default = 1
% 
% Outputs:
%   yes: 1 = item is between them, 0 = not
% 
% John Grogan, 2020

%% 

if ~exist('isInclusive','var') || isempty(isInclusive)
    isInclusive = 1; 
end


%%

if isInclusive
    
    yes = data >= limits(1) & data <= limits(2);
    
else
    
    yes = data > limits(1) & data < limits(2);
    
end


end