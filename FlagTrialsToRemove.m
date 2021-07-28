function [EEG, com] = FlagTrialsToRemove(EEG, varargin)
% insert a flag into epochs that occur within trials marked as toRemove for
% artefact rejection purposes.
%
% Inputs:
%   EEG: EEGLAB data structure, containing EEG.result.data.toRemove
%   Flag     : value of flag(s) for (0 to 16) to set as trigger.
%               (0 = remove flag). default = 1
%   Channel  : index of channel to insert flag into, default = last one.
%
% Outputs:
%   EEG      : EEGLab data structure, with flags inserted for trials with
%   toRemove
%
% Example:
%   EEG = FlagTrialsToRemove(EEG, 'Flag', [5 6], 'Channel', 64);
%
% Requires:
%   EEGLAB toolbox (https://sccn.ucsd.edu/eeglab/index.php)
%   ERPLab toolbox (https://erpinfo.org/erplab)
%   parsepvpairs from Matlib (https://github.com/sgmanohar/matlib)
%   nancat from Matlib
%
% Written by John P. Grogan, 2020.
%

%% error check inputs

com = ''; % setup for fails
if nargin<1
    error('no inputs. see "help FlagTrialsToRemove"');
end
if isobject(EEG) % eegobj
    whenEEGisanObject % calls a script for showing an error window
    return
end

%% if there are multiple datasets

if length(EEG) > 1
    [ EEG, com ] = eeg_eval( 'FlagTrialsToRemove', EEG, 'warning', 'on');
    return;
end

%% inputs

if ~isfield(EEG, 'result') || ~isfield(EEG.result.data(1), 'toRemove')
    error('EEG must have EEG.result.data.toRemove');
end

toRemove = [EEG.result.data.toRemove];

if length(toRemove) ~= EEG.trials
    error('toRemove does not match number of trials in EEG');
end


params = {'Flag', 'Channel'}; % param names
values = {1, EEG.nbchan}; % default values

[flag, chanArray] = parsepvpairs(params, values, varargin{:});

if all( flag(:) ~= [0:16],'all' )
    error('Flag must be matrix of 0-16. see "help FlagTrialsToRemove"');
end

%% for each trigger/flag



% %% find all triggers inside window in each epoch
% 
% % nancat will give lots of warnings about empty cells, which can be ignored
% % here, so I am disabling that warning - requies nancat from Feb 2020
% warning('OFF', 'NANCAT:emptyCells');
% 
% % get all events into matrix
% allEvent = nancat(1, EEG.epoch.eventtype); % cell array [nEpochs x nEvents]
% allEventLat = nancat( nancat(1, EEG.epoch.eventlatency) ); % get latencies too
% 
% warning('ON', 'NANCAT:emptyCells'); % turn this warning back on
% 
% 
% % replace empty cells with empty chars, so regexp will work
% allEvent(cellfun(@isempty, allEvent)) = {''};
% 
% % find the triggers - usign regexpi on cell
% eventFlags = ~cellfun(@isempty, regexpi(allEvent,triggers{j})); % set the flag for each
% 
% % find events within window
% withinWindow = allEventLat >= tWindow(1) & allEventLat <= tWindow(2);
% 
% % set those to zero
% eventFlags(~withinWindow) = 0;


%% insert flags for each epoch into EEG structure

EEG = flagEpoch(EEG, find(toRemove), flag, chanArray);

%% print number of rejected

perreject = mean(toRemove) * 100; % how many we rejected
fprintf('FlagTrialsToRemove() rejected: %.1f%% of total trials %s (flag = %d).\n', perreject, flag);





%% display artefact table

pop_summary_AR_eeg_detection(EEG, ''); % show table at the command window
EEG = eeg_checkset( EEG );


%% return the command?

com = 'EEG = FlagTrialsToRemove( EEG'; % start of command

pv = varargin; % copy
for i = 1:length(varargin)
    if mod(i,2) % param
        com = [com sprintf(', ''%s''', num2str(varargin{i}))];
        
    else % value
        if ischar(varargin{i})
            if ~strcmpi(varargin{i},'off')
                com = [com sprintf( ', %s', varargin{i})];
            end
        else
            if iscell(varargin{i})
                if isstr(varargin{i}{1})
                    fn2resstr = vect2colon(cell2mat(cellfun(@str2num, varargin{i},'UniformOutput',0)), 'Sort', 'On');
                    fnformat = '%s';
                else
                    fn2resstr = vect2colon(cell2mat(varargin{i}), 'Sort','on');
                    fnformat = '{%s}';
                end
            else
                fn2resstr = mat2colon(varargin{i}, 'Sort','on');
                fnformat = '%s';
            end
            com = [com sprintf( [', %s ' '%s' fnformat], fn2resstr)];
        end
    end
end
com = [com '); % ' datestr(now) ]; % print command and timestamp

% append to history
history = cellstr(EEG.history);
EEG.history = char([history; {com}]);

EEG = eeg_checkset( EEG ); % check again, will also format history

end

