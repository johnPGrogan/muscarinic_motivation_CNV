function EEG = DeletePracticeTrialTriggers(EEG)
% function EEG = DeletePracticeTrialTriggers(EEG)
% deletes the triggers for the practice trial(s) (if exist)
%
% Should be run after merging files, but before saccade triggers are
% inserted, and before DeleteTriggersForTrials()
%

if ~isfield(EEG, 'result')
    EEG = InsertMatFileIntoEEG(EEG, behDataFolder);
end
result = EEG.result;

% make a copy of EEG for editing
EEG1 = EEG;

%% check for practice trials

if isfield(result, 'practiceResult') && ~isempty(result.practiceResult)
    nPracTrials = length(result.practiceResult); % number of practice trials
else
    warning('no practice trials were found, so cannot be deleted');
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

if isfield(result, 'deleted') && any(result.deleted)
    nDelTrials = length(result.deleted);
else
    nDelTrials = 0;
end

% check numbers of trials
if nTrigTrials == nDataTrials
    fprintf('number of trigger trials matches data trials, not deleting any\n');
    return;
    
else
    fprintf('different number of trials from result.data and triggers\n');
    
    if nTrigTrials == (nDataTrials + nPracTrials)
        fprintf('the number of extra trials matches the number of practice trials, so will be deleted\n');
        
    elseif nTrigTrials  == (nDataTrials + nPracTrials + nDelTrials)
        fprintf('the number of extra trials match number of practice + number to be deleted. deleting practice trials only for now\n');
        
    elseif nTrigTrials  == (nDataTrials + nDelTrials)
        fprintf('the number of extra trials matches with the number in result.deleted. not deleting any here, these will be deleted in DeleteTriggersForTrials\n');
        return;
        
    end
    
end

%% delete practice trials

indsToDel = []; % indices of triggers to remove

% i will assume practice trials are at the start only

for i = 1:nPracTrials
    
    ev = itiInds(i); % event to start deleting from
    
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
    if newNTrigTrials ~= (nDataTrials + nDelTrials)
        warning('new number of trigger trials does not match with number of data + toDelete trials');
        keyboard;
    end
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