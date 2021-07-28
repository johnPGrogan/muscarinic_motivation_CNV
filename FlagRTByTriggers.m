function EEG = FlagRTByTriggers(EEG, varargin)
% function EEG = FlagRTByTriggers(result, varargin)
% flag trials where RT (from triggers) is outside given limits. Will find
% RT in trial, even if outside the epoch
% 
% Inputs:
%   EEG = eeglab structure
%   varargin: param-values of things to filter by:
%       'rtLimits'   = [min max] RT limits (inclusive) to filter by. You 
%           can use [0 Inf] to not filter at all.
%       'targRespTriggers' = 1x2 cell of regexpi strings for finding
%           target & response triggers. default = {'48|51', '54'};     
%       'flag'       = flag(s) to insert (e.g. [1 9]). default is [1]
% 
% Outputs:
%   EEG
% 

% default inputs
params = {'rtLimits', 'targRespTriggers', 'flag'};
values = {[0 Inf],{'48|51', '54'}, [1]};

% get inputs
[rtLimits, targRespTriggers, flag] = parsepvpairs(params, values, varargin{:});

if isempty(rtLimits); fprintf('rtLimts is empty, will not flag by RT\n'); end



%% flag triggers in continuous data

% trigString = '27|30|45|72|75'; % eye stuff or pause

% look in rest of trial, outside of epoch

% get times of epoch event
% look between that and start/end of trial

% flag any triggers within that period

trialLimTrigs = {'^6$', '^6$|66|75|78|84'}; % start of trial, start of next | end of this trial/block/expt

% get event indices of bin events
eventInds = arrayfun(@(x) x.event(find([x.eventbini{:}] > 0,1,'first')), EEG.epoch);

% get urevent
% EEG.urevent(EEG.event(eventi(1)).item)

events = {EEG.urevent.type}; % original events, before epoching 

limits(:,1) = cellRegexpi(events, trialLimTrigs{1}); % index of each itiStart
limits(:,2) = cellRegexpi(events, trialLimTrigs{1}); % index of each trial end (or start of next trial)

nEvents = length(eventInds);

%% flag RT

[urInd, startInds, endInds, targInds, respInds] = deal(zeros(nEvents,1));
trialTrigs = cell(nEvents,1);

% get target/dist trigger for this trial/epoch
for i = 1:nEvents
    
    urInd(i,1) = EEG.event(eventInds(i)).item; % index of this urevent
    
    % find urevent index of start:end
    startInds(i,1)  = find(limits(1:urInd(i),1),1,'last');
    endInds(i,1)    = min([find(limits(urInd(i):end,2),1,'first') + urInd(i) - 1, length(events)]);
    
    trialTrigs{i,1} = {EEG.urevent(startInds(i):endInds(i)).type}; % all triggers in this trial
  
    % find targ/resp triggers on trial
    targInds(i,1) = find(cellRegexpi(trialTrigs{i}, targRespTriggers{1})) + startInds(i) - 1; % target
    respInds(i,1) = find(cellRegexpi(trialTrigs{i}, targRespTriggers{2})) + startInds(i) - 1; % response fixation
    
end

% get latencies of each
targLat = [EEG.urevent(targInds).latency]';
respLat = [EEG.urevent(respInds).latency]';

% take difference - convert into ms
rt = (respLat - targLat) ./EEG.srate * 1000;

% compare to limits
toFlag = rt<rtLimits(1) | rt>rtLimits(2);

% insert those flags
EEG = flagEpoch(EEG, find(toFlag), flag, 1);


end
