function validInds = FindTriggersBetweenTriggers(trigs, dim, trigsToFind, windowTrigs, indsToReturn, strict, emptyVals)
% function validInds = FindTriggersBetweenTriggers(trigs, trigsToFind, windowTrigs, dim, indsToReturn, [strict], [emptyVals])
% 
% find triggers between two other triggers. if before/after not given,
% uses the beginning/end of epoch. if criteria is 0 it will allow triggers
% that are before 'before' but where 'after' does not occur within epoch,
% or similar (default is 1).
% 
% Inputs:
%   trigs = cell array of strings containing triggers to search within
%   dim = dimension to search along
%   trigsToFind = string or cell string containing a single trigger to
%       find using regexpi
%   windowTrigs = cell array of 2 trigger strings, between which the
%       trigger must lie. if one is empty string ('') and strict==1, it 
%       will use the beginning/end of epoch, as appropriate, otherwise it
%       will find nothing.
%   indsToReturn = 'all' = return all triggers within window, 'first' =
%       first trigger only, 'last' = last trigger only, N (scalar) = Nth
%       trigger only. Any epochs with fewer than N valid triggers will
%       return NaN.
%   strict = 1 excludes trials where either windowTrig is not present, 0 =
%       allows these trials, using the beginning/end of epoch. default = 1
%   emptyVals = value to use for when no valid indices are found (default =
%       [], can also use NaN or 0).
% 
% Outputs:
%   validInds = cell array containing valid indices for each trial,
%       matching the size of trigs. Epochs without valid triggers will have
%       NaN instead
% 
% John Grogan, 2020.

%% inputs checking
if ~(iscell(trigs) && iscellstr(trigs))
    error('trigs must be a cell array of strings');
end

if isempty(dim) || dim>ndims(trigs) || mod(dim,1)~=0
    error('dim must be an integer <= ndims(trigs)');
end

if ~((iscell(trigsToFind) && iscellstr(trigsToFind) && ~cellfun(@isempty, trigsToFind)) || (ischar(trigsToFind) && ~isempty(trigsToFind)))
    error('trigsToFind must be a string or cell string');
end

if ~(iscell(windowTrigs) && length(windowTrigs)==2 && iscellstr(windowTrigs))
    error('windowTrigs must be a 2 cell array of strings');
end

if ~((ischar(indsToReturn) && any(strcmp(indsToReturn, {'first','last','all'}))) || (isnumeric(indsToReturn) && length(indsToReturn)==1))
    error('indsToReturn must be a string ("first", "last", or "all") or a scalar index');
end

if ~exist('strict','var') || isempty(strict)
    strict = 1;
end

if ~exist('emptyVals','var')
    emptyVals = [];
end

%% find the indices of triggers within trials

trigInds = apply(dim, @(x) {find(x)}, cellRegexpi(trigs, trigsToFind)); % find indices of trigger within trial
trigInds(cellfun(@isempty, trigInds)) = {0}; % replace empties

startWindowInds = apply(dim, @(x) {find(x)}, cellRegexpi(trigs, windowTrigs{1})); % find indices of trigger within trial
endWindowInds = apply(dim, @(x) {find(x)}, cellRegexpi(trigs, windowTrigs{2})); % find indices of trigger within trial

% some windows can have multiple triggers, so just use first?
startWindowInds(cellfun(@length, startWindowInds)>1) = cellfun(@(x) x(1), startWindowInds(cellfun(@length, startWindowInds)>1),'UniformOutput',0);

%% does trigger occur between before+after?

if strict % strict = before + after must be present
    
    % replace empties so that NO triggers fall within them
    startWindowInds(cellfun(@isempty, startWindowInds)) = {Inf}; % 
    endWindowInds(cellfun(@isempty, endWindowInds)) = {0}; % replace

    validInds = cellfun(@(x,y,z) x(x>y & x<z), trigInds, startWindowInds, endWindowInds,'UniformOutput',0); 

    
else % lax = before or after can be absent
    
    % replace empties so that ALL triggers fall within them
    startWindowInds(cellfun(@isempty, startWindowInds)) = {0}; % 
    endWindowInds(cellfun(@isempty, endWindowInds)) = {Inf}; % replace

    validInds = cellfun(@(x,y,z) x(x>y & x<z), trigInds, startWindowInds, endWindowInds,'UniformOutput',0); 

end


%% return number of inds requested

switch indsToReturn
    case 'first'
        validInds(~cellfun(@isempty, validInds)) = cellfun(@(x) x(1), validInds(~cellfun(@isempty, validInds)),'UniformOutput',0);
    case 'last'
        validInds(~cellfun(@isempty, validInds)) = cellfun(@(x) x(end), validInds(~cellfun(@isempty, validInds)),'UniformOutput',0);
    case 'all'
        % return all inds
    otherwise % number
        validInds(cellfun(@length, validInds) < indsToReturn) = {emptyVals}; % 
        validInds(cellfun(@length, validInds) < indsToReturn) = cellfun(@(x) x(indsToReturn), validInds(cellfun(@length, validInds) < indsToReturn),'UniformOutput',0);
end


%% replace empty

validInds(cellfun(@isempty, validInds)) = {emptyVals};
end

