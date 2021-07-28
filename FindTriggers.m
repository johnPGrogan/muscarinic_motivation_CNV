function trigInds = FindTriggers(channelData, threshold)
% trigInds = FindTriggers(channelData, threshold)
% finds indices of where the channelData changes by more than threshold
% (positive and negative change)
% and returns those indices
% Inputs:   channelData is a vector of the voltages from the EEG channel
%               which had the trigger inputs
%           threshold is the voltage change threshold to use for finding
%               the triggers
% Outputs:  trigInds is a matrix with columns for indices of onset and
%               offset of each trigger detected
% 
% John Grogan, 2019.


trigOn = diff(channelData) > threshold; % find where V changes (onset)
trigOff = diff(channelData) < -threshold; % and goes back (offset)

trigOnInds = findregions(trigOn); % find all those changes
trigOffInds = findregions(trigOff); % and for offsets

trigInds = [trigOnInds(:,1), trigOffInds(:,1)]; % return indices of onset and offsets

end