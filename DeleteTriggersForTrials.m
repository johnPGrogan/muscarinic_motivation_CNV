function EEG = DeleteTriggersForTrials(EEG)
% function EEG = DeleteTriggersForTrials(EEG)
% deletes the triggers for trials corresponding to result.deleted (which
% were removed during readAllEDF due to restarting expt etc.
%
% Should be run after merging files, but before saccade triggers are
% inserted
%

if ~isfield(EEG, 'result')
    EEG = InsertMatFileIntoEEG(EEG, behDataFolder);
end
result = EEG.result;

% make a copy of EEG for editing
EEG1 = EEG;

%% are there result.deleted?

if isfield(result, 'deleted') && any(result.deleted)
    nDelTrials = length(result.deleted);
else
    warning('there are no trials to be deleted in result.deleted');
    return;
end

%% get stims

% get any boundary events (from merging two datasets)

cellregexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

isText = ~cellregexpi({EEG1.event.type}, '\d')'; % find stim without numbers in them

if any(isText) % pull these out, and replace later
    
    textEvents = EEG1.event(isText);
    EEG1.event(isText) = [];
    
end

% numeric stimuli
stims(:,1) = (nancat(1, EEG1.event.latency) - 1) ./ EEG1.srate; % latency into time
stims(:,2) = str2num(nancat(1,EEG1.event.type)); % trigger code
stims(:,3) = ones(size(stims,1),1); % set to 1,


%% get number of stim trials

trigs = stims(:,2); % get trigs

itiInds = find(trigs==6); % index of ITI (start of trial)

nTrigTrials = length(itiInds);
nDataTrials = length(result.data);

if isfield(result, 'practiceResult') % were there practice trials?
    nPracTrials = length(result.practiceResult);
else
    nPracTrials = 0;
end

% check they add up
if nTrigTrials ~= nDataTrials
    fprintf('different number of trials from result.data and triggers\n');
    
    if nTrigTrials == (nDataTrials + nDelTrials)
        fprintf('the number of extra trials match with result.deleted, so will be deleted\n');
        
    elseif nTrigTrials == (nDataTrials + nDelTrials + nPracTrials) || nTrigTrials == (nDataTrials + nPracTrials)
        warning('there are still practice trials in here - please run DeletePracticeTrialTriggers\n');
        keyboard;
    end
    
end

indsToDel = []; % indices of triggers to remove

% for each trial to delete
for i = 1:length(result.deleted)
    tr = result.deleted(i); % trial in triggers
    ev = itiInds(tr); % event to start deleting from
    
    % find when next trial/block/expt starts/ends
    j = find(any(trigs(ev+1:end) == [3 6 78 81 84],2), 1, 'first');
    %             for j = 1:length(itiInds)
    %                 if any(trigs(ev+j) == [3 6 78 81 84])
    %                     break;
    %                 end
    %             end
    
    % store indices to delete
    indsToDel = [indsToDel ev:(ev+j - 1)];
end

% disp and wait
disp('deleting:')
[a,b] = CountUnique(trigs(indsToDel));[b,a]

% delete them
trigs(indsToDel,:) = [];

% check it
disp('leaving:')
[a,b] = CountUnique(trigs);[b,a]
CheckTrigCodeNumbers(trigs); % check different triggers match up 
keyboard;

stims(indsToDel,:) = []; % actually remove from stims


%% check the new trial numbers add up

newNTrigTrials = sum(stims(:,2)==6); % new number of ITIs

if newNTrigTrials ~= nDataTrials
    warning('new number of trigger trials does not match with number of data trials');
    keyboard;
end


%% put these stims back into EEG
% make sure still in order
stims = sortrows(stims); % order by time

EEG1.event = struct('type',cellfun(@num2str,num2cell(stims(:,2)),'UniformOutput',0), 'latency',num2cell(stims(:,1)*EEG1.srate+1), 'urevent',num2cell([1:length(stims)]'), 'duration', num2cell(zeros(length(stims),1)) );
EEG1.urevent = struct('type',cellfun(@num2str,num2cell(stims(:,2)),'UniformOutput',0), 'latency',num2cell(stims(:,1)*EEG1.srate+1));
EEG1 = eeg_checkset( EEG1 );


%% if there were any text events

if any(isText)
    
    EEG1.event = [EEG1.event; textEvents]; % append
    
    % sort back into time order
    [~,i] = sort([EEG1.event.latency]);
    EEG1.event = EEG1.event(i);
    
end

% return EEG1 as EEG if it made it this far
EEG = EEG1;