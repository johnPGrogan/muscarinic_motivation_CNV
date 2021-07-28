function EEG = InsertMatFileIntoEEG(EEG, behDataFolder, mFileName)
% insert mat file name into EEG
% and the matfile
% Inputs:
%   EEG = eeglab structure
%   behDataFolder = path where mFile is
%   mFileName [optional] = name of mat file to load, will be taken from EEG
%   name info if not given


%% inputs

if ~exist('mFileName','var') || isempty(mFileName)
    if isfield(EEG, 'matName') && ~isempty(EEG.matName)
        mFileName = EEG.matName;
    else
        mFileName = EEG.setname(5:32); % get mat name
    end
end

if ~strcmp(mFileName(end-3:end), '.mat') %% append extension if not already
    mFileName = [mFileName '.mat'];
end

%% load mat file

mFile = load(fullfile(behDataFolder, mFileName)); % load mat file


%% store name in EEG

% check if already there
if isfield(EEG, 'matName')
    if ~equals(EEG.matName, mFileName)
        fprintf('EEG already contains a different matName (%s) to the new one (%s)\n', EEG.matName, mFileName);
        keyboard;
    end 
end

EEG.matName = mFileName; % store name
    
%% store result structure

% check if already there
if isfield(EEG, 'result')
    if ~equals(EEG.result, mFile.result)
        fprintf('EEG already contains a different result to the new one \n');
        disp(EEG.result);
        disp(mFile.result);
        keyboard;
    end 
end


EEG.result = mFile.result; % store data


end
    
    
    
