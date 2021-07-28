function [EEG, com] = FlagTriggersInWindow(EEG, varargin)
% Flag any triggers occuring within a time window of each epoch, for
% artefact rejection purposes.
%
% Inputs:
%   EEG: EEGLAB data structure
%   Triggers : trigger values to flag - cell of strings, or vector of numbers
%   TWindow  : vector of time window (ms) [start end]. default is entire epoch
%   Flag     : value of flag for (0 to 16) to set for each trigger. One row
%                   of flags per trigger. (0 = remove flag). default = 1
%   Channel  : index of channel to insert flag into, default = last one.
%
% Outputs:
%   EEG      : EEGLab data structure, with flags inserted
%
% Example:
%   EEG = FlagTriggersInWindow(EEG, 'Triggers', [31 32], 'TWindow', [-1500 0], 'Flag', [5 6], 'Channel', 64);
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
    error('no inputs. see "help FlagTriggersInWindow"');
end
if isobject(EEG) % eegobj
    whenEEGisanObject % calls a script for showing an error window
    return
end


%% parse arguments

% define defaults
params = {'Triggers', 'TWindow', 'Flag', 'Channel'}; % param names
values = {NaN, NaN, 1, EEG.nbchan}; % default values

% assign param-values into variables - function from matlib package
[triggers, tWindow, flag, chanArray] = parsepvpairs(params, values, varargin{:});

% check for errors
if ~iscell(triggers) && all(isnan(triggers))
    error('Triggers parameter not defined. see "help FlagTriggersInWindow"');
    
elseif all(isnumeric(triggers)) % is it a vector of numbers?
    triggers = cellstr(num2str(triggers'))'; % make into cell of strings
else %if not vector of numbers
    if ~iscell(triggers) % if not a cell array
        error('Triggers are not numbers or string of chars. see "help FlagTriggersInWindow"');
    else
        x = cell2mat(triggers'); % make matrix
        if all(isnumeric(x)) % is it a cell array of numbers?
            triggers = cellstr(num2str(triggers'))'; % make into cell of strings
        elseif all(ischar(x)) && isnumeric(str2double(x)) % is it a cell array of string numbers?
            % do nothing
        else
            error('Triggers are not numbers or string of chars. see "help FlagTriggersInWindow"');
        end
    end
end

if any(size(tWindow) ~= [1 2])
    error('TWindow must be a two-element row vector, e.g. [-1000 0]. see "help FlagTriggersInWindow"');
    
elseif tWindow(2) <= tWindow(1)
    error('TWindow ends before it starts - elements must increase. see "help FlagTriggersInWindow"');
    
elseif all( flag(:) ~= [0:16],'all') || any(size(flag,1) ~= size(triggers,2))
    error('Flag must be matrix of 0-16, with 1 row per Trigger. see "help FlagTriggersInWindow"');
    
end


%% if there are multiple datasets

if length(EEG) > 1
    options1 = reshape([params; [{triggers}, {tWindow}, {flag}, {chanArray}] ], 1, []); % make into 1 row
    [ EEG, com ] = eeg_eval( 'FlagTriggersInWindow', EEG, 'warning', 'on', 'params', options1);
    return;
end


%% for each trigger/flag

for j = 1:length(triggers)
    %% find all triggers inside window in each epoch

    % nancat will give lots of warnings about empty cells, which can be ignored
    % here, so I am disabling that warning - requies nancat from Feb 2020
    warning('OFF', 'NANCAT:emptyCells');

    % get all events into matrix
    allEvent = nancat(1, EEG.epoch.eventtype); % cell array [nEpochs x nEvents]
    allEventLat = nancat( nancat(1, EEG.epoch.eventlatency) ); % get latencies too

    warning('ON', 'NANCAT:emptyCells'); % turn this warning back on


    % replace empty cells with empty chars, so regexp will work
    allEvent(cellfun(@isempty, allEvent)) = {''};

    % find the triggers - usign regexpi on cell
    eventFlags = ~cellfun(@isempty, regexpi(allEvent,triggers{j})); % set the flag for each

    % find events within window
    withinWindow = allEventLat >= tWindow(1) & allEventLat <= tWindow(2);

    % set those to zero
    eventFlags(~withinWindow) = 0;


    %% insert flags for each epoch into EEG structure

    flaggedEpochs = find(any(eventFlags,2)); % which epochs have any
    isRT = 0; % maybe allow this to change

    for i = 1:length(flaggedEpochs)

        % put a flag into this EEG.epoch. can insert all channels flags at once
        [EEG, errorm]= markartifacts(EEG, flag(j,:), chanArray, 1:length(chanArray), flaggedEpochs(i), isRT);
        if errorm==1
            error(['ERPLAB says: There was not latency at the epoch ' num2str(i)])
        elseif errorm==2
            error('ERPLAB says: invalid flag (0<=flag<=16)')
        end

    end

    % if we didn't reach the last one, do it now, so that EEG.reject is correct size
    if length(EEG.reject.rejmanual) ~= EEG.trials
        
        EEG.reject.rejmanual(end+1:EEG.trials) = 0; % fill with zeros      
        EEG.reject.rejmanualE(:, end+1:EEG.trials) = 0; 
        
    end
    %% print number of rejected

    if j==1
        fprintf('FlagTriggersInWindow() rejected: \n');
    end
    perreject = mean(any(eventFlags,2)) * 100; % how many we rejected
    fprintf('%.1f%% of total trials for trigger %s (flag = %d).\n', perreject, triggers{j}, flag(j));

end



%% display artefact table

pop_summary_AR_eeg_detection(EEG, ''); % show table at the command window
EEG = eeg_checkset( EEG );


%% return the command?

com = 'EEG = FlagTriggersInWindow( EEG'; % start of command

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

