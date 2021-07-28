function EEG = myEyeArtefactFlagging(EEG, varargin)
% function EEG = myEyeArtefactFlagging(EEG, varargin)
% wrapper function to call my eye artefact flagging routine, uses eeglab.
% Flags muscle artefacts using pop_artextval, and uses FlagTriggersInWindow
% to remove blinks/saccades (detected by eyetracker and trigger flagged) in
% the same window.
% appends '_artflag' to EEG.setname
% Inputs:
%   EEG = eeglab structure
%   varargin: pairs of 'param', value
%       artWindow = time window for artefact flagging in epoch
%       save = 1 = save
%       saveFolder = folder to save in
% 
% Output: EEG structure with flags set, and also stored in
% EEG.reject.rejmanual
% 

%% parse inputs

params = {'artWindow', 'save', 'saveFolder', 'Threshold'};
values = {[], 0, [], [ -200 200]};

[artWindow, isSave, saveFolder, thresh] = parsepvpairs(params, values, varargin{:});

if isempty(artWindow)
    disp('no artWindow given, so will use entire epoch')
    windargs = {};
else
    windargs = {'TWindow', artWindow};
end


nChannels = EEG.nbchan; % assumes last 4 channels are hEOG, vEOG, elX, elY

%% muscle - voltage

EEG  = pop_artextval( EEG , 'Channel',  1:nChannels-2, 'Flag',  [1 2],...
     'Threshold', thresh, windargs{:} ); % GUI: 19-Aug-2019 12:57:11


%% EOG
% 
% % blinks - moveing window - EOG up/down (VEOG)
% EEG  = pop_artmwppth( EEG , 'Channel',  63, 'Flag',  [1 3],...
%     'Threshold',  70, 'Twindow', artWindow, 'Windowsize',  200, 'Windowstep',...
%     100 ); % GUI: 19-Aug-2019 13:04:16
% 
% 
% % eye movements - step-like - HEOG
% EEG  = pop_artstep( EEG , 'Channel',  62,  'Flag',  [1 4],...
%     'Threshold',  25, 'Twindow', artWindow, 'Windowsize',  200,...
%     'Windowstep', 50 ); % GUI: 19-Aug-2019 13:14:47
% 

%% also run on eye tracking synched data
% % > 1vis deg or diameter of fix cross
% 
% % blinks - moveing window - EOG up/down (VEOG)
% EEG  = pop_artmwppth( EEG , 'Channel',  nChannels-1:nChannels,...
%     'Flag',  [1 5], 'Threshold',  200, 'Twindow', artWindow,...
%     'Windowsize',  200, 'Windowstep', 100 ); % GUI: 19-Aug-2019 13:04:16
% 
% 
% % eye movements - step-like - HEOG
% EEG  = pop_artstep( EEG , 'Channel',  nChannels-1:nChannels,...
%     'Flag',  [1 6], 'Threshold', 50, 'Twindow', artWindow,...
%     'Windowsize',  100, 'Windowstep', 30 ); % GUI: 19-Aug-2019 13:14:47
% 


%% exclude eyetracker detected blinks or saccades within artWindow

% insert flags for blink start/end, saccade start/end
EEG = FlagTriggersInWindow(EEG, 'Triggers', [31 32 52 53],...
    'twindow', artWindow, 'Flag', [1 7; 1 8; 1 8; 1 8]);


%% flag channels and epochs by variance

% ignore eye-tracking traces
%   chanVar = nanvar(EEG(iPP,iSes).data(1:end-2,:,:), 1, [2 3]); % variance within and between epochs

%   epochVar = sq(nanvar(EEG(iPP,iSes).data(1:end-2,:,:), 1, [1 2])); % variance within epochs, across channels

%% pull together for later rejection

EEG = eeg_rejsuperpose( EEG, 1, 1, 0, 0, 0, 0, 0, 1); % manual only


% reset all
% EEG  = pop_resetrej( EEG , 'ArtifactFlag',  1:16, 'ResetArtifactFields', 'on', 'UserFlag',  1:16 )
%% check & save

EEG.setname = [EEG.setname '_artflag'];
EEG = eeg_checkset( EEG );

if isSave
    if isempty(saveFolder); error('saveFolder not given');end
    
    EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',saveFolder);
end


end