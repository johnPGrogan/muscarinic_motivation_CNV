function ok = CheckEEGMatData(EEG, kb, varargin)
% function ok = CheckEEGMatData(EEG, kb, varargin)
% checks the dates of the EEG file, mat file, eyelink files all match
% checks these against the records
% will give errors/warnings if there are mismatches
% Input:
%       EEG = eeglab structure, with EEG.result = result
%       kb = use keyboard instead of error when mismatches (will pause parfor)
%       varargin: strings of which checks to call: 'ppID', 'dates', 'RT', 'trialtype', 'triggerRTs'
% 
% Calls subfunctions for CheckMatEdfDates, CheckMatSaccTimes, 
% CheckTriggerTrialTypes, CheckTriggerRTs
% 
% Output: ok = vector. 1 = returned, 0 = errors, NaN = not run.
% 
%% is there a mat file?

if ~isfield(EEG, 'matName') && ~isfield(EEG, 'result')
    error('no matName or result found in EEG');
end
if ~exist('kb','var') || isempty(kb)
    kb = 1;
end

ok = NaN(1, length(varargin));



%% run specified functions

if any(strcmpi(varargin, 'ppID'))
    ok(strcmpi(varargin, 'ppID')) = CheckEEGPpInfo(EEG, kb); % check ppID and sesNo
end

if any(strcmpi(varargin, 'dates'))
    ok(strcmpi(varargin, 'dates')) = CheckMatEdfDates(EEG.result, kb); % check date of mat file against edf
end

if any(strcmpi(varargin, 'RT'))
    ok(strcmpi(varargin, 'RT')) = CheckMatSaccTimes(EEG.result, kb); % compares RT in mat & edf
end

if any(strcmpi(varargin, 'trialTypes'))
    if isfield(EEG,'epoch') && ~isempty(EEG.epoch) % epoch needs different code
        ok(strcmpi(varargin, 'trialTypes')) = CheckTriggerTrialTypesEpoched(EEG, kb); % compares trigger and mat trial type info
    else
        ok(strcmpi(varargin, 'trialTypes')) = CheckTriggerTrialTypes(EEG, kb); % compares trigger and mat trial type info
    end
end

if any(strcmpi(varargin, 'triggerRT'))
    if isfield(EEG,'epoch') && ~isempty(EEG.epoch)
        ok(strcmpi(varargin, 'triggerRT')) = CheckTriggerRTsEpoched(EEG, kb); % compares trigger and mat RTs
    else
        ok(strcmpi(varargin, 'triggerRT')) = CheckTriggerRTs(EEG, kb); % compares trigger and mat RTs
    end
end
    

end


function ok = CheckEEGPpInfo(EEG, kb)
%% check ppID and sesNo

    ok = 0;

    i = regexpi(EEG.setname, 'TRI_\d\d\d');
    eegID = EEG.setname(i:i+6);
    eegSess = EEG.setname(i+8);

    if ~strcmpi(eegID, EEG.result.params.ppID)
        KbOrError(kb, 'eegID does not match EEG.result.params.ppID', {eegID, EEG.result.params.ppID});
    end

    if str2num(eegSess) ~= EEG.result.params.sessionNumber
        KbOrError(kb, 'eeg session number  does not match EEG.result.params.sessionNumber', {eegSess, EEG.result.params.sessionNumber});
    end
    
    ok = 1;

end



function ok = CheckMatEdfDates(result, kb)
% check edf file(s) dates match mat file dates
% Input: result structure from RunExperiment
% 
% will give warnings/errors if there are mismatches

ok = 0;
%% check number of files
nEdf = length(result.edfFiles); % number of edf files
nMat = length(result.startTimes);

if nEdf ~= nMat
    KbOrError(kb, 'different number of edf and mat files. nEdf = %d, nMat = %d', [nEdf, nMat]);
end

%% parse the startTimes date to check that too

for i = 1:nMat
    matDates{i} = result.startTimes{i}([6 7 9 10 12 13 15 16]);
end

%% check that each edf file matches matDates

for i = 1:nEdf
    if ~equals(result.edfFiles{i}(1:8), matDates{i})
        disp(result);
        fprintf('result.edfFiles dates do not match result start times\n')
        
        if (str2double(matDates{i}) - str2double(result.edfFiles{i}(1:8))) == 1
            fprintf('timestamps are out by 1 min: %s, %s\n', result.edfFiles{i}, matDates{i});
            
        elseif equals(result.edfFiles{i}(1:4), matDates{i}(1:4))
            fprintf('the dates do match: %s, %s\n', result.edfFiles{i}, matDates{i});
            
        else
            KbOrError(kb, 'the dates do not match either', {result.edfFiles{i}, matDates{i}});
        end
    end
end

       
ok = 1;
end


function ok = CheckMatSaccTimes(result, kb)
% check RTs in result.s against result.data
% if they differ by >10ms, it will pause and show the differences
% (actually takes time from ITI to Response, as the codes are the same so
% it is simpler, rather than comparing RT)

ok = 0;
%% is there a result.s?

if ~isfield(result,'s') || isempty(result.s)
    KbOrError(kb,'there is no result.s',{});
end


%% do the dimensions match up?

if any(size(result.s) ~= size(result.data))
    KbOrError(kb, 'result.s and result.data have different numbers of trials',[]);
end

%% compare the times from ITI to Response

% milliseconds
saccTimes = ([result.s.saccadeaccepted_t]-[result.s.startITI_t]); % from eye tracker
behTimes = ([result.data.saccadeaccepted]-[result.data.startITI]) .*1000; % from matlab

[maxDiff, maxInd] = max(abs(saccTimes - behTimes)); % biggest difference?

if sum(abs(saccTimes-behTimes)>10)>2 % milliseconds
    

    figure();
    subplot(1,2,1)
    plot(behTimes, saccTimes,'x');
    xlabel('behRT (ms)');
    ylabel('saccRT (ms)');

    subplot(1,2,2);
    plot( saccTimes - behTimes, 'x' )
    ylabel('saccRT - behRT (ms)');
    
    
    KbOrError(kb, 'saccRT and behRT do not match', {});
end
    
ok = 1;
end

function ok = CheckTriggerTrialTypes(EEG, kb)
% check trigger code trial types against result.data
% Input: EEG (eeglab structure, containing EEG.result = result)
% 
% It will give warnings and errors if things do not match up
% 

ok = 0;
%% do we have the info?

if ~isfield(EEG, 'event') || ~isfield(EEG, 'result')
    KbOrError(kb,'EEG must contain event and result',{});
end

result = EEG.result; % extract for simplicity

%% get the trigger info

events = {EEG.event.type}; % get all event codes

iti = strcmp(events, '6'); % ITI codes

if sum(iti) ~= size(result.data, 2)
    KbOrError(kb, 'number of ITI codes does not match number of trials', {sum(iti), size(result.data,2)});
end

% trial type info to compare
trigCodes = struct( 'targetDir',    {{'9','12'; -1, 1}},...
                    'fixLoc',       {{'15', '18', '21'; 1 2 3}},...
                    'reward',       {{'36', '39'; 0 50}},...
                    'distractor',   {{'48', '51'; 0 1}} );

fn = fieldnames(trigCodes); % trigs to find

for i = 1:length(fn)

    tmp = NaN(size(events,2), length(trigCodes.(fn{i}))); % init

    for j = 1:length(trigCodes.(fn{i})) % is there any trigger code?

        % put the new value in (same as from result.data)
        tmp( strcmp(events, trigCodes.(fn{i}){1,j}),j) = trigCodes.(fn{i}){2,j};
    end

    trig.(fn{i}) = nansum( tmp( any(~isnan(tmp),2),:),2 )'; % get only those trigs, in order presented

end


%% get the trial types from beh data

behNames = {'targetDir', 'previousTargetPos','rewardMax','distractor'};

for i = 1:length(behNames)
    beh.(fn{i}) = [result.data.(behNames{i})];
end


%% compare

for i = 1:length(fn)

    if any( trig.(fn{i}) ~= beh.(fn{i}) )
        KbOrError(kb, 'trigger codes do not match trial types', {fn{i}, [trig.(fn{i});  beh.(fn{i})]});
    end
end


ok = 1;
    

end


function ok = CheckTriggerTrialTypesEpoched(EEG, kb)
% check trigger code trial types against result.data if EEG is epoched
% Input: EEG (eeglab structure, containing EEG.result = result)
% 
% It will give warnings and errors if things do not match up
% 

ok = 0;
%% do we have the info?

if ~isfield(EEG, 'event') || ~isfield(EEG, 'result')
    KbOrError(kb, 'EEG must contain event and result', {});
end

result = EEG.result; % extract for simplicity

%%
if length(EEG.epoch) ~= size(result.data,2)
    KbOrError(kb, 'number of epochs does not match number of trials', {length(EEG.epoch), size(result.data,2)});
end

% check if event codes match

% trial type info to compare
trigCodes = struct( 'targetDir',    {{'^9','12'; -1, 1}},...
                    'fixLoc',       {{'15', '18', '21'; 1 2 3}},...
                    'reward',       {{'36', '39'; 0 50}},...
                    'distractor',   {{'48', '51'; 0 1}} );

fn = fieldnames(trigCodes); % trigs to find

events = nancat(1,EEG.epoch.eventtype);
events(cellfun(@isempty, events)) = {'a'}; % replace empties so regexpi works
for i = 1:length(fn)

    tmp = NaN(size(events,1), length(trigCodes.(fn{i}))); % init

    for j = 1:length(trigCodes.(fn{i})) % is there any trigger code?

        % put the new value in (same as from result.data)
        tmp( any(cellRegexpi(events, trigCodes.(fn{i}){1,j})~=0,2), j) = trigCodes.(fn{i}){2,j};
    end

    trig.(fn{i}) = nansum( tmp( any(~isnan(tmp),2),:),2 )'; % get only those trigs, in order presented

end

% get beh trigs
behNames = {'targetDir', 'previousTargetPos','rewardMax','distractor'};
for i = 1:length(behNames)
    beh.(fn{i}) = [result.data.(behNames{i})];
end


% compare non-empties
for i = 1:length(fn)

    if ~isempty(trig.(fn{i})) && length(trig.(fn{i}))==length(beh.(fn{i})) && any( trig.(fn{i}) ~= beh.(fn{i}) )
        KbOrError(kb, 'trigger codes do not match trial types', {fn{i}, [trig.(fn{i});  beh.(fn{i})]});
    end
end


ok = 1;

end


function ok = CheckTriggerRTs(EEG, kb)
% check trigger code RTs match with behavioural RTs
% Input: EEG (eeglab structure, containing EEG.result = result)
% 
% It will give warnings and errors if things do not match up
% 

ok = 0;
%% do we have the info?

if ~isfield(EEG, 'event') || ~isfield(EEG, 'result')
    KbOrError(kb, 'EEG must contain event and result',{});
end

result = EEG.result; % for ease

%% get the trigger info

events = {EEG.event.type}; % get all event codes

stim = find(strcmp(events, '48') | strcmp(events, '51')); % ITI codes

if size(stim,2) ~= size(result.data, 2)
    KbOrError(kb, 'number of stim codes does not match number of trials', {size(stim,2), size(result.data,2)});
end


% response
resp = find(strcmp(events, '54'));

if size(resp,2) ~= size(stim,2)
    KbOrError(kb,'number of response codes does not match number stim codes', {size(resp,2), size(stim,2)});
end


% get times
respTimes = [EEG.event(resp).latency] ./EEG.srate; % convert into seconds
stimTimes = [EEG.event(stim).latency] ./EEG.srate;

eegRT = (respTimes - stimTimes) *1000;

%% get the RT from behavioural data

behRespTimes = [result.data.saccadeaccepted] ; % from matlab, seconds
behStimTimes = [result.data.starttarget] ;

behRT = (behRespTimes - behStimTimes) * 1000 ; 


%% compare

[maxDiff, maxInd] = max(abs(behRT - eegRT)); % biggest difference?

if sum(abs(behRT-eegRT)>10)>2 % ms
    
    figure();
    subplot(1,2,1)
    plot(behRT, eegRT,'x');
    xlabel('behRT (ms)');
    ylabel('eegRT (ms)');

    subplot(1,2,2);
    plot( behRT - eegRT,'x')
    ylabel('behRT - eegRT (ms)');
    
    
    KbOrError(kb, 'behTimes and eegTimes do not match on trial', maxInd);

end
    

ok = 1;
end

function ok = CheckTriggerRTsEpoched(EEG, kb)
% check trigger code RTs match with behavioural RTs in epoched data
% Input: EEG (eeglab structure, containing EEG.result = result)
% 
% It will give warnings and errors if things do not match up
% 

ok = 0;
%% do we have the info?

if ~isfield(EEG, 'epoch') || ~isfield(EEG, 'result')
    KbOrError(kb, 'EEG must contain epoch and result',[]);
end

result = EEG.result; % for ease

%% get the trigger info

events = nancat(1, EEG.epoch.eventtype); % get triggers per epoch
events(cellfun(@isempty, events)) = {'a'}; % replace empties so regexpi works

if size(events,1) ~= size(result.data, 2)
    KbOrError(kb, 'number of epochs does not match number of trials', {size(events,1), size(result.data,2)});
end


stimInds = cellRegexpi(events, '48|51') ~= 0;
respInds = cellRegexpi(events, '54') ~= 0;
eventTimes = nancat(1, EEG.epoch.eventlatency);

[stimTimes, respTimes] = deal(NaN(size(events,1),1));
for i = 1:size(stimInds,1)
    stimInd = find(stimInds(i,:),1);
    respInd = find(respInds(i,:),1);
    
    if ~isempty(stimInd)
        stimTimes(i,1) = eventTimes{i,stimInd};
    end
    if ~isempty(respInd)
        respTimes(i,1) = eventTimes{i,respInd};
    end
end

% get times
% respTimes = respTimes; % convert into seconds
% stimTimes = stimTimes;

eegRT = (respTimes - stimTimes);

%% get the RT from behavioural data

behRespTimes = [result.data.saccadeaccepted] ; % from matlab, seconds
behStimTimes = [result.data.starttarget] ;

behRT = (behRespTimes - behStimTimes)' * 1000 ; 


nMissing = sum(isnan(eegRT) & behRT<max(EEG.times)); % missing eeg rt which should be within epoch
if nMissing > 0
    KbOrError(kb, 'there are %d trials with missing eeg RTs, likely because RT was greater than epoch', nMissing);
end


%% compare

[maxDiff, maxInd] = max(abs(behRT - eegRT)); % biggest difference?

if sum(abs(behRT-eegRT)>10)>2 % seconds
    
    
%%
    figure();
    subplot(1,2,1)
    plot(behRT, eegRT,'x');
    xlabel('behRT (ms)');
    ylabel('eegRT (ms)');

    subplot(1,2,2);
    plot( behRT - eegRT,'x')
    ylabel('behRT - eegRT (ms)');
  %%  
    
    KbOrError(kb, 'behTimes and eegTimes do not match', maxInd);
end
    

ok = 1;
end

function KbOrError(kb, msg, disps)
% if kb, use warning+keyboard, else use error
if kb
    warning(msg);
    disp(disps);
    keyboard;
else
    disp(disps);
    error(msg);
end

end