function EEG = trioEEGProcessLoad(ID, sectionsToRun)
% trioEEGProcessLoad(ID, sectionsToRun)
% Run this on each person's session individually. 
% sectionsToRun determines which parts are run (e.g. can omit convert/load
% etc to save time if already have been done).
% Inputs: 
%   ID = string with ppID and session number e.g. 'Tri_001_1'
%        If there are two files for a session, use a string that is
%        common to both.
%   sectionsToRun = struct of bools on which sections to run, must contain
%       these fields, with bool (0 or 1) for each:
%           'convert',  1,... % convert ov2mat
%           'load',     1,... % trioEEGLoad (incl trigger stuff)
%           'merge',    1,... % whether to remerge split files
%           'delTrigs', 1,... % delete trigs from practice & deleted trials
%           'saccTrig', 1,... % insert sacc/blink trigs
%           'check',    1,... % check info, times, trials match
% 
% Saves files during convert, load, merge. Saves EEG file at end, and
% returns EEG structure also
% 

if ~exist('ID','var') || ~isstr(ID)
    error('ID must be a string file name ID');
end

if ~exist('sectionsToRun','var') || ~all(isfield(sectionsToRun, {'convert','load','merge','delTrigs','saccTrig','check'}))
    error('sectionsToRun not supplied correctly. see help trioEEGProcessLoad');
end

addpath ../../Experiments/MyFuncs/

%% check eeglab is open, open if not

hh = findobj('tag', 'EEGLAB');
if isempty(hh); eeglab; end % call this to setup path and initialise

%% set up folders

triFolder = 'D:/groganj/Work/Trihexy/';
if ~exist(triFolder,'file')
    triFolder = 'D:/Backups/Trihexy/';
end

eegDataFolder = [triFolder 'EEGData/']; % where to load EEG data from
behDataFolder = [triFolder 'BehData/']; % where EDF and Mat files are 
analysisFolder = [triFolder 'EEGAnalysis/']; % where to save analysis files

%% get file names

f = dir(eegDataFolder); % find files in data folder 
fileNames = {f.name}';

eegFileInds = ~cellfun(@isempty, regexpi(fileNames, ['eeg_TrioRewardCue3_EEG_' ID])); % get the ones that match ID
eegFiles = fileNames(eegFileInds); 

fprintf('%d eeg files found:\n',length(eegFiles)) % show
disp(eegFiles)

% remove the extension
eegFiles = cellfun(@(x) x(1:regexp(x,'\.')-1), eegFiles, 'UniformOutput',0);
eegFiles = unique(eegFiles); % remove duplicates (in case .ov and .mat present)

% do the same for behavioural mat file 
f = what(behDataFolder);
fileNames = f.mat;
mFileInds = ~cellfun(@isempty, regexpi(fileNames, ['TrioRewardCue3_EEG_' ID]));
mFiles = fileNames(mFileInds);


fprintf('%d result mat files found:\n',length(mFiles)) % show
disp(mFiles)

%% check we have any files

if isempty(eegFiles) || isempty(mFiles)
    error('eegFiles or mFiles is empty. see above');
end

%% convert the .ov to .mat

for i = 1:length(eegFiles) % each eeg file separately

    if sectionsToRun.convert || ~exist(fullfile(eegDataFolder, [eegFiles{i} '.mat']), 'file') && ~exist(fullfile(analysisFolder, [eegFiles{i} '.set']), 'file')
        convert_ov2mat(fullfile(eegDataFolder, [eegFiles{i} '.ov']), fullfile(eegDataFolder, [eegFiles{i} '.mat']));
    end
    
end

%% get edf file name for eye-EEG trigger parsing and synching later

for i = 1:length(mFiles) % if multiple mat files give
  mFile = load(fullfile(behDataFolder, mFiles{i})); % load it

  elFileName1 = mFile.result.edfFiles; % get edf file names
  
  if ~iscell(elFileName1); elFileName1 = {elFileName1};end % put into cell
  
  for j = 1:length(elFileName1) % get each edf file name
    elFileNames{j,i} = elFileName1{j}(1:end-4);
  end
  
end

elFileNames = reshape(elFileNames,[],1); %columnise

fprintf('%d edf files found:\n',length(elFileNames)) % print filenames
disp(elFileNames)

if any(size(elFileNames) ~= size(eegFiles)) % if they don't match up
    disp(elFileNames)
    disp(eegFiles)
    disp('different number of eyelink .edf files and eeg .mat files found. Please correct this or abort')
    keyboard; % 
end


%% load .mat into eeg

for i = 1:length(eegFiles) % for each eeg file
    
    if sectionsToRun.load || ~exist(fullfile(analysisFolder, [eegFiles{i} '.set']), 'file') || ~exist(fullfile(analysisFolder, [eegFiles{i} '.fdt']), 'file')
        trioEEGLoad( fullfile(eegDataFolder, [eegFiles{i} '.mat']), analysisFolder, eegFiles{i}, behDataFolder, elFileNames{i}, 3, mFiles{1});
    end
    
end


%% merge two datasets - if experiment was stopped partway through

if length(eegFiles)>1 % if there are two files

    eegNames = nancat(1,eegFiles{:}); % get their names
    commonChars = eegNames(1, all(diff(eegNames,[],1)==0,1)); % pick common characters
    
    % merge
    if sectionsToRun.merge || ~exist(fullfile(analysisFolder, [commonChars(1:end-1) '.set']), 'file') || ~exist(fullfile(analysisFolder, [commonChars(1:end-1) '.fdt']), 'file')
        fprintf('merging two datasets\n');
        
        for i = 1:length(eegFiles) % load each file
            ALLEEG(i) = pop_loadset('filename',[eegFiles{i} '.set'],'filepath',analysisFolder);
            ALLEEG(i) = eeg_checkset( ALLEEG(i) );
        end

        EEG = pop_mergeset( ALLEEG, [1  2], 0); % merge them - set up for two files for now
        EEG = eeg_checkset( EEG );


        EEG.setname = commonChars(1:end-1); % set name based on common characters
        EEG.filename = EEG.setname;
        
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',analysisFolder); % save
    
    end
end
    

%% load up file if not already present

if ~exist('EEG','var') || ~isfield(EEG,'data') || isempty(EEG.data)
    if exist('commonChars','var')
        EEG = pop_loadset('filename',[commonChars(1:end-1) '.set'], 'filepath', analysisFolder);
    else
        EEG = pop_loadset('filename',[eegFiles{1} '.set'], 'filepath', analysisFolder);
    end
end


%% load matlab file with data too

EEG = InsertMatFileIntoEEG(EEG, behDataFolder, mFiles{1});


%% delete trials that were removed

if sectionsToRun.delTrigs
    EEG = DeletePracticeTrialTriggers(EEG); % remove practice trials
    
    EEG = DeleteTriggersForTrials(EEG); % remove missing/deleted trials
end

%% insert saccades & blinks, if not already present
if sectionsToRun.saccTrig
    
    EEG = InsertSaccadeTriggersIntoEEG(EEG, behDataFolder);
    
end

%% check data all matches up
if sectionsToRun.check
   
    % checks ppID, sesNo, dates of edf/mat files, rts, trial types, trigRTs,
    % all match up
    ok = CheckEEGMatData(EEG, 'ppID', 'dates', 'RT', 'trialTypes', 'triggerRT');
    if ~all(ok==1)
        warning('not all checks were run');
        disp(ok);
        keyboard;
    end

end

%% update name and save

EEG.setname = [EEG.setname '_processed'];

EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',analysisFolder);

end