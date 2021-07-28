function [postSaccTimes, saccLats] = GetTimesBeforeTriggers(saccCodes, beforeCode, trigLats, trigs, times)
% function [postSaccTimes, saccLats] = GetTimesBeforeTriggers(saccCodes, beforeCode, trigLats, trigs, times)
% Returns boolean matrix of whether each time sample was after the saccade
% trigger codes given, on each trial.
% Requires: FindTriggersBetweenTrials.m
%
% Inputs: 
%   saccCodes = cell array of strings of trigger codes to find
%   beforeCode = string of trigger before which all triggers are excluded
%       (e.g. beforeCode='B' will exclude any triggers occuring before the
%       binning code)
%   trigLats = [nTr,nPP,nSes,nTrigs] matrix of trigger latencies
%   trigs = [nTr,nPP,nSes,nTrigs] cell matrix of trigger strings
%   times = [1 nTimes] vector of times (EEG.times)
% 
% Outputs:
%   postSaccTimes = [1 nTimes nTr nPP nSes nSaccCodes] boolean matrix of
%       whether each time sample is after each saccadeCode within each
%       trial.
%   saccLats = [nTrials, nPPs, nSes, nSaccCodes] cell array of latencies
% 
% John Grogan, 2020

% saccCodes = {'52','53','54'}; % 52=start, 53=end, 54=accepted

%% input checking
if ~( iscell(saccCodes) && iscellstr(saccCodes) )
    error('saccCodes must be a cell array of regexp strings');
end

if ~( (iscell(beforeCode) && iscellstr(beforeCode)) || (isstr(beforeCode)))
    error('beforeCode must be a string or cell string trigger');
end

if ~( isnumeric(trigLats))
    error('trigLats must be numeric matrix');
end

if ~( iscell(trigs) && all(size(trigLats) == size(trigs)) && iscellstr(trigs) )
    error('trigs must be a cell array of trigger strings, with '''' replacing missing values');
end

if ~( isnumeric(times) )
    error('times must be a [1 x N] vector of times');
end

%%
nSaccCodes = length(saccCodes); % number of codes to find

[nTrials, nPP, nSes, nTrigs] = size(trigs); % sizes

% convert trigLats into [nTr, nPp, nSes] cell mat
trigLatMat = mat2cell(trigLats,ones(nTrials,1),ones(1,nPP),ones(1,nSes),nTrigs); % sq
trigLatMat = cellfun(@sq, trigLatMat,'UniformOutput',0); % squeeze

% % find target/distr index
% targVals = cellRegexpi(trigs, beforeCode);
% targInds = apply(4, @(x) {find(x)}, targVals); % find saccInds within trial
% % targLats = cellfun(@(x,y) x(y), trigLatMat, targInds,'UniformOutput',0); % get latencies
% targInds(cellfun(@isempty, targInds)) = {0}; % replace empty
% targInds = cellfun(@(x) {x(1)}, targInds); % get first targ

% saccInds = NaN(nTrials, nPP, nSes, nSaccCodes);    
saccLats = cell(nTrials,nPP,nSes,nSaccCodes);
for i = 1:nSaccCodes % get saccade indices

    % find indices of saccades between beforeCode and end of trial, replacing missing with []
    saccIndsi = FindTriggersBetweenTriggers(trigs, 4, saccCodes{i}, [beforeCode {''}], 'first', 0, []);
    
%     % find saccCodes indices
%     saccVals = cellRegexpi(trigs, saccCodes{i});
%     saccIndsi = apply(4, @(x) {find(x)}, saccVals); % find saccInds within trial
%     % get first one after target
%     saccIndsi = cellfun(@(x,y) x(find(x>y,1,'first')), saccIndsi, targInds,'UniformOutput',0); 
    
    saccLats(:,:,:,i) = cellfun(@(x,y) x(y), trigLatMat, saccIndsi,'UniformOutput',0); % latencies
%     saccIndsi(cellfun(@isempty, saccIndsi)) = {NaN}; % replace empty with NaN
%     saccInds(:,:,:,i) = nancat(saccIndsi); % un-cell
    
end

% get indices to NaN (at the time of the triggers)
postSaccTimes = cellfun(@(x) times>min([x Inf]), saccLats,'UniformOutput',0); % 1 = after
% [1 nTimes nTr nPP nSes nSaccCodes]