function EEG = myEpoching(EEG, varargin)
% EEG = myEpoching(EEG, varargin)
% wrapper function for epoching data. calls pop_creabasiceventlist,
% pop_binlister, pop_epochbin. It can either baseline to a period within
% the epoch defined (using eeglab), or can baseline to a period within a
% separate epoch (using separate bdf file), or it can baseline to a
% period defined by a trigger earlier in the current epoch (manually).
% 
% inputs:
%   EEG = eeglab data structure
%   optional param-value pairs:
%       resFolder = folder to save eventlist and file in
%       bdfFile   = bin descriptor file for epoching e.g. './trio_EEG_targ_BDF.txt'
%       timeWindow = time window to epoch e.g. [-1500 0]
%       baselineMethod = string saying which type of baselining to use.
%           (default is 'window')
%               'none' = no baselining
%               'window' = baseline to a window relative to epoch trigger.
%                   Must provide baselineWindow
%               'withinEpoch' = baseline to a window set by different
%                   trigger within the same epoch. Must provide
%                   baselineTriggers, baselineWindow and timeWindow2
%               'separateEpoch' = baseline to a window in a separate set of
%                   epochs, defined by a new bdf file. Must provide
%                   baselineBDFFile and baselineWindow
%       baselineWindow = baseline window. if using baselineMethod of
%           'window', it will be applied relative to the event trigger in
%           the bdfFile. If using 'withinEpoch', it will be applied
%           relative to the baselineTriggers. If using 'separateEpoch', it
%           will be relative to the event trigger set by baselineBdfFile
%       baselineBDFFile = BDF filename to use for epoching, if different to
%           main epoching file
%       baselineTriggers = regexp string of triggers to use for defining a 
%           baseline period within the epochs set by bdfFile
%       timeWindow2 = smaller window than timeWindow to use after doing a
%           manual baselining within epochs (using
%           baselineTriggersInEpoch), relative to the trigger event in
%           bdfFile.
%       suffix = suffix to append to EEG.setname
%       save = 1 = save file (optional)
%
% Outputs:
%   EEG structure epoched
%   optionally saves the EEG structure also

%% param-values

params = {'resFolder', 'bdfFile', 'timeWindow', 'baselineMethod', 'baselineWindow', 'baselineBdfFile', 'baselineTriggers',  'timeWindow2', 'suffix', 'save'};
values = {[], [], [], 'window', [], [], [], [], [], 0};

[resFolder, bdfFile, timeWindow, baselineMethod, baselineWindow, baselineBdfFile,  baselineTriggers,  timeWindow2, suffix, isSave] = parsepvpairs(params, values, varargin{:});

if isempty(resFolder);  error('resFolder must be specified');   end
if isempty(bdfFile);    error('bdfFile must be specified');     end
if isempty(timeWindow); error('timeWindow must be specified');  end

% input checking on baseline method
if isempty(baselineMethod)% [] supplied
    error('baselineMethod is empty');
elseif strcmpi(baselineMethod, 'none') && (~isempty(baselineWindow) && ~strcmpi(baselineWindow, 'none')) % says 'none' but a baselineWindow is given different to 'none'
    error('baselineMethod is none, but a baselineWindow is supplied');
elseif strcmpi(baselineMethod, 'window') && ~isempty(baselineBdfFile) % says 'window' but extra bdf file
    error('baselineMethod is window, but a baselineBdfFile is supplied');
elseif strcmpi(baselineMethod, 'withinEpoch') && ~isempty(baselineBdfFile) % says withinEpoch but baselineBdfFile
    error('baselineMethod is withinEpoch, but a baselineBdfFile is supplied');
elseif strcmpi(baselineMethod, 'withinEpoch') && (isempty(baselineTriggers) || isempty(baselineWindow) || isempty(timeWindow2))
    error('baselineMethod is withinEpoch, but must supply baselineTriggers, baselineWindow and timeWindow2');
elseif strcmpi(baselineMethod, 'separateEpoch') && (isempty(baselineBdfFile) || isempty(baselineWindow))
    error('baselineMethod is separateEpoch, must supply baselineBdfFile and baselineWindow');
end


%% CREATING EPOCHS - main epoch

% create event list
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', [resFolder EEG.setname '_eventlist.txt'] ); % GUI: 27-Jun-2019 09:28:35
EEG = eeg_checkset( EEG );
% EEG.setname = [EEG.setname '_elist'];

% create bins using the BDF file
EEG  = pop_binlister( EEG , 'BDF', bdfFile, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' ); % GUI: 27-Jun-2019 09:29:33
EEG = eeg_checkset( EEG );
% EEG.setname = [EEG.setname '_bins'];

%% do baseline manually?

% if there is no separate bdf file
switch baselineMethod 
    
    case 'none' % no baseline, just epoch
        
        fprintf('no baselining\n');
        
        % epoch data and apply baseline
        EEG = pop_epochbin( EEG , timeWindow, 'none'); % GUI: 27-Jun-2019 09:29:57
        
        EEG.baseline = [];
        
    case 'window' % baseline to window relative to event trigger
    
        fprintf('baselining within main file\n');

        % epoch data, no baseline for now
        EEG = pop_epochbin( EEG , timeWindow,  'none'); % GUI: 27-Jun-2019 09:29:57
    
        % manually get baseline
        EEG.baseline = nanmean(EEG.data(:,isBetween(EEG.times, baselineWindow),:),2);
        
        EEG.data = EEG.data - EEG.baseline; % subtract
        
    case 'withinEpoch' % baseline relative to different trigger in same epoch

        fprintf('manually baselining within epoch with window [%i %i] relative to trigger %s\n', baselineWindow(1), baselineWindow(2), baselineTriggers);
        
        % epoch data. no baseline
        EEG = pop_epochbin( EEG , timeWindow,  'none'); % GUI: 27-Jun-2019 09:29:57
        
        % manually baseline to earlier trigger within these epochs, trim
        % epochs. store the baseline in EEG.baseline
        EEG = baselineWithinEpoch(EEG, baselineTriggers, baselineWindow, timeWindow2); % do so 
    
        
    case 'separateEpoch' % baseline to triggers in new epochs
    
        fprintf('manually baselining from separate epochs using bdf:%s, with window [%i %i]\n', baselineBdfFile, baselineWindow(1), baselineWindow(2));

        % recursion - new epochs, no baseline
        EEGBaseline = myEpoching(EEG, 'resFolder', resFolder, 'bdfFile', baselineBdfFile, 'timeWindow', baselineWindow, 'baselineMethod', 'none');

        % get main epoch
        EEG = pop_epochbin( EEG , timeWindow,  'none'); % GUI: 27-Jun-2019 09:29:57

        % check sizes match up
        if any([EEG.nbchan EEG.trials size(EEG.data,1) size(EEG.data,3)] ~= [EEGBaseline.nbchan EEGBaseline.trials size(EEGBaseline.data,1) size(EEGBaseline.data,3)])
            warning('EEG and EEGBaseline epoch data sizes do not match up: %s', EEG.setname);
%             keyboard;
        end


        % get mean within baseline period
        preCueBaseline = nanmean(EEGBaseline.data, 2);

        % subtract
        EEG.data = EEG.data - preCueBaseline;
        
        % store baseline
        EEG.baseline = preCueBaseline;
    
    
    otherwise
    
        error('error with baselineBdfFile');
end


%% check

EEG = eeg_checkset( EEG );

EEG.setname = [EEG.setname suffix];


%% save

if isSave

    EEG = pop_saveset( EEG, 'filename', EEG.setname,'filepath',resFolder);
    
end

end


function EEG = baselineWithinEpoch(EEG, baselineTriggers, baselineWindow, timeWindow2)
% function EEG = baselineWithinEpoch(EEG, baselineTriggers, baselineWindow, timeWindow2)
% epoch to large window, manually find baseline period within this,
% baseline with that, then re-epoch to trim window. Store baseline in
% EEG.baseline
  
    % does every epoch contain the bsTrigger?
    tc = nancat(1, EEG.epoch.eventtype);
    bsTrig = zeros(size(tc));
    for i = 1:numel(tc)
        if ~isempty(tc{i})
            bsTrig(i) = cellRegexpi(tc(i),baselineTriggers)>0;
        end
    end
    if ~all(sum(bsTrig,2)==1)
        warning('not all epochs have the baseline trigger: %s',EEG.setname);
%         keyboard;
    end
    
    % find the baseline period within that epoch
    bsLat = nancat(1, EEG.epoch.eventlatency); % get latencies
    bsLat(~bsTrig) = {0}; % only keep bsTrig
    bsLat = nancat(bsLat); % matrix
    
    bsTrigLat = sum(bsLat, 2); % sum to keep bsTrigLat
    
    % get the activity in baseline
    baseline = NaN(EEG.nbchan,1, EEG.trials);
    for i = 1:EEG.trials
        % get indices of baseline window
        [p1, p2] = window2sample(EEG, bsTrigLat(i) + baselineWindow, EEG.srate, 'relaxed');
        
        % get activity within
        baseline(:,1,i) = nanmean(EEG.data(:,p1:p2,i),2);
    end
    
    % subtract baseline
    EEG.data = EEG.data - baseline;
    
    % store 
    EEG.baseline = baseline;
    
    % trim that epoch - by new time
    EEG = pop_epochbin( EEG, timeWindow2,  'none');
    
end