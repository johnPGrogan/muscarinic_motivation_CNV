function EEG = FlagTriggersInTrial(EEG, varargin)
% function EEG = FlagTriggersInTrial(result, varargin)
% flag trials for removal based on triggers in EEG, works across entire
% trial, not just the epoch (use FlagTriggersInWindow for within-epoch
% flagging).
% The trigger flagging looks at a period defined by two triggers (e.g.
% trial start + end). The RT flagging works on targ/resp triggers.
%
% Inputs:
%   EEG = eeglab structure
%   varargin: param-values of things to filter by:
%       'triggers'   = regexpi string of triggers to exlude from entire
%           trial.
%       'trialLimTrigs' = 1x2 cell array, each containing a regexpi string
%           used to define the start/end triggers for trial (or section).
%       'flag'       = flag(s) to insert (e.g. [1 9]). default is [1]
%
%
% Outputs:
%   EEG
%

% default inputs
params = {'triggers', 'trialLimTrigs', 'flag'};
values = {'27|30|45|72|75', {'^6$', '^6$|66|72|78|84'}, 1};

% get inputs
[trigString, trialLimTrigs, flag] = parsepvpairs(params, values, varargin{:});

if isempty(trigString); error('triggers is empty'); end
if isempty(trialLimTrigs); error('trialLimTrigs is empty'); end


%% flag triggers in continuous data
% find triggers in continuous data correspond to start/end of this trial
% flag any of the trigString that occur within that period

% trigString = '27|30|45|72|75'; % eye stuff or pause
%
% start of trial, start of next | end of this trial/block/expt
% trialLimTrigs = {'^6$', '^6$|66|75|78|84'};

% get event indices of bin events
eventInds = arrayfun(@(x) x.event(find([x.eventbini{:}] > 0,1,'first')), EEG.epoch);%, 'UniformOutput',0);

% get urevent
events = {EEG.urevent.type}; % original events, before epoching

% indices of start/end
limits(:,1) = cellRegexpi(events, trialLimTrigs{1}); % index of each itiStart
limits(:,2) = cellRegexpi(events, trialLimTrigs{1}); % index of each trial end (or start of next trial)

%% find all events in trial

nEvents = length(eventInds);
[urInd, startInds, endInds, toFlag] = deal(zeros(nEvents,1));
trialTrigs = cell(nEvents,1);
for i = 1:nEvents
    
    urInd(i,1) = EEG.event(eventInds(i)).item; % index of this urevent
    
    % find urevent index of start:end
    startInds(i,1)  = find(limits(1:urInd(i),1),1,'last'); % last 'start' before this epoch
    endInds(i,1)    = min([find(limits(urInd(i):end,2),1,'first') + urInd(i) - 1, length(events)]); % first 'end' after this epoch
    
    trialTrigs{i,1} = {EEG.urevent(startInds(i):endInds(i)).type}; % all triggers in this trial
    
    toFlag(i,1) = any(cellRegexpi(trialTrigs{i}, trigString)); % any to flag?
    
end

% insert those flags
EEG = flagEpoch(EEG, find(toFlag), flag, 1);





end


