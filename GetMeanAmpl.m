function ampl = GetMeanAmpl(data, times, timeWindow)
% function ampl = GetMeanAmpl(data, times, timeWindow)
% 
% Return mean amplitude within time windows
% 
% Inputs:
%   data: matrix with time in 2nd dimension, any other shape allowed
%   times: [1 nTimes] vector of times
%   timeWindow: [min max] times of window (inclusive)
% 
% Outputs:
%   ampl: matrix of amplitudes, matching size of data but with size(ampl,2)
%       as 1
% 
% John Grogan, 2020

%% check inputs

if all(size(timeWindow) ~= [1 2]) || min(timeWindow) < min(times) || max(timeWindow) > max(times)
    error('timeWindow must be a [1 x 2] vector, within the limits of times variable');
end

if size(times,2) ~= size(data,2)
    error('times must be the same size as size(data,2)');
end

%% get window

timeInds = isBetween(times, timeWindow, 1); % inclusive

% set all times outside window to NaN, then average over dim2
% this allows us to average within window regardless of size of data
sizes = size(data);
sizes(2) = 1;
timeInds = repmat(timeInds, sizes);

data(~timeInds) = NaN; 

%% get mean

ampl = nanmean(data, 2);