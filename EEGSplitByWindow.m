function eegByWindow = EEGSplitByWindow(data, times, windows)
% function eegByWindow = EEGSplitByWindow(data, times, windows)
% Split the EEG data by time windows, for taking SD across time within a
% window
% Inputs:
%   data = EEG data [chans, time, bin, pp, ses, trial]
%   times = [1 N] vector of times
%   windows = vector of window boundaries
% 
% Outputs:
%   eegByWindow = matrix of eeg data split by window [pp, window, bin, sess, chan, trial, sample]
% 
% John Grogan, 2020

%% window SD
nSamples = length(times); % number of samples in epoch


nWindows = length(windows) - 1;
windowInds = NaN(1,nSamples);
for i = 1:nWindows
    windowInds(times >= windows(i) & times<windows(i+1)) = i;
end

% get SD across window within each trial
eegByWindow = permute(groupMeans(data,2,windowInds,'dim'),[7,6,2,3,4,5,1]); % [pp, window, bin, sess, chan, trial, sample]

end
