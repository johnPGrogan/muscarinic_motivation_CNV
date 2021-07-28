function EEG = InsertSaccadeTriggersIntoEEG(EEG, behDataFolder, minSaccSize, blinksToo)
% call InsertSaccadesAsTriggers without having to do entire trioEEGLoad
% script. it will insert them right into the EEG structure
% 

if ~exist('minSaccSize','var') || isempty(minSaccSize)
    minSaccSize = 1 / .0203; % size of min saccade to accept in pixels
end
if ~exist('blinksToo', 'var') || isempty(blinksToo)
    blinksToo = 1; % include blink triggers too
end

%% get any boundary events (from merging two datasets)

cellregexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

isText = ~cellregexpi({EEG.event.type}, '\d')'; % find stim without numbers in them

if any(isText) % pull these out, and replace later
    
    textEvents = EEG.event(isText);
    EEG.event(isText) = [];
    
end 

%% extract stims from EEG events

stims(:,1) = (nancat(1, EEG.event.latency) - 1) ./ EEG.srate; % latency into time
stims(:,2) = str2num(nancat(1,EEG.event.type)); % trigger code
stims(:,3) = ones(size(stims,1),1);

if any(stims(:,2)==[31 32 52 53],'all')
    warning('there are already blink or saccade triggers inserted into this EEG event structure. duplicates will be removed');
end

%% load up result.s

if ~isfield(EEG, 'result')
    EEG = InsertMatFileIntoEEG(EEG, behDataFolder);
end
result = EEG.result;

%% check readAllEDF is stored in result.s

if length(result.data) ~= length(result.s) || ~isfield(result, 's') 
    edf2ascMat( fullfile(behDataFolder, result.edfFiles) );
    result.s = readAllEDF(result, behDataFolder);
    disp('you may want to save this')
    keyboard;
%     save( fullfile(behDataFolder, mFile), 'result', '-append');
end
%% get triggers for saccades & blinks

stims2 = InsertSaccadesAsTriggers(result, stims, minSaccSize, blinksToo); % include blinks

%% insert these back in

EEG.event = struct('type',cellfun(@num2str,num2cell(stims2(:,2)),'UniformOutput',0),'latency',num2cell(stims2(:,1)*EEG.srate+1),'urevent',num2cell([1:length(stims2)]'),'duration', num2cell(zeros(length(stims2),1)));
EEG.urevent = struct('type',cellfun(@num2str,num2cell(stims2(:,2)),'UniformOutput',0),'latency',num2cell(stims2(:,1)*EEG.srate+1));
EEG = eeg_checkset( EEG );


%% if there were any text events

if any(isText)
    
    EEG.event = [EEG.event; textEvents]; % append
    
    % sort back into time order
    [~,i] = sort([EEG.event.latency]);
    EEG.event = EEG.event(i);
    
end

%% save?


end