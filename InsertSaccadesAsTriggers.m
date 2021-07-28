function stimsNew = InsertSaccadesAsTriggers(result, stimsOld, minSaccSize, blinksToo)
% InsertSaccadesAsTriggers(result, stimsOld, minSaccSize, blinksToo)
% detects saccades and blinks from the edf file for each trial, inserts the
% starts and ends as triggers, based on time since ITI. Also inserts a
% trigger for the first saccade after the target onset (or the saccade
% acceptance time if no saccades detected after target onset).
% Inputs:
%   result = behavioural data, including result.s = readAllEdf()
%   stimsOld = triggers from InsertTriggerCodes or similar, [time, code, duration]
%   minSaccSize = minimum saccade size to keep (in pixels)
%   blinksToo = whether to do blinks also [1 = default];
% 
% Outputs:
%   stimsNew = combined output with stimsOld and saccade onset/offset
%   triggers included, in ascending time order
% 
% John Grogan, 2020

%%  snip saccades from start ITI to end trial

% turn off this warning
warning('off','NANCAT:emptyCells');

% don't interpolate blinks - will find blinks from NaNs
[raw, info]=snipSaccades(result.s,'startITI_t','endfeedback_t', 'minsize', minSaccSize, 'interpolate',0);
% there may be 1 fewer as practice are not snipped but are in stims2

% turn warning back on
warning on;
%% get times of starts of saccades from start of recording

itiTimes = stimsOld( stimsOld(:,2)==6,1); % time ITI trigger was sent, in seconds
if length(result.s) ~= length(itiTimes)
    if isfield(result, 'practiceResult') && length(itiTimes) == (length(result.s) + length(result.practiceResult))
        itiTimes(1) = []; % remove practice
    else
        disp('itiTimes and result.s do not match up. manually adjust these and type dbcont to continue');
        keyboard;
    end
end

saccStartTimes = itiTimes + info.sRT/1000; % get start times of saccades in seconds from start of recording

% check trials don't overlap
finalSaccStart=alignRight(saccStartTimes); % find last saccade in each trial
if any(finalSaccStart(1:end-1,end) > itiTimes(2:end))
    disp('final saccade starts later than beginnning of next trial')
    keyboard;
end

%% get times of ends of saccades 

saccEndTimes = saccStartTimes + info.sDur/1000; % end of saccades

% check trials don't overlap
finalSaccEnd=alignRight(saccEndTimes); % find last saccade in each trial
if any(finalSaccEnd(1:end-1,end) > itiTimes(2:end))
    disp('final saccade end later than beginnning of next trial')
    keyboard;
end


%% blinks
if ~exist('blinksToo','var') || blinksToo % if doing blinks
    
    % find blinks from non-interpolated traces
    blinkInds = apply(2, @(x) {findregions(x)}, isnan(raw));
    blinkInds = permute(nancat(blinkInds),[3,2,1]);

    % remove trailing nans
    blinkInds(repmat(blinkInds(:,2,:) >= size(raw,2),1,2)) = NaN;

    blinkStartTimes = sq(blinkInds(:,1,:)); % take start indices (ms from ITI trigger)
    blinkEndTimes = sq(blinkInds(:,2,:)); % end indices

    % check that there are the same number of each
    if any(isnan(blinkStartTimes) & ~isnan(blinkEndTimes), 'all')
        disp('blink starts and ends do not match up')
        keyboard;
    end

    blinkDurations = blinkEndTimes - blinkStartTimes; % get duration (ms)

    % maybe threshold these blink durations? 
    minBlinkDur = 100; % ms

    validBlinks = blinkDurations >= minBlinkDur; % only keep ones that reach this
    blinkStartTimes(~validBlinks) = NaN; % exclude others
    blinkEndTimes(~validBlinks) = NaN;
    blinkDurations(~validBlinks) = NaN;

    % convert to seconds, and add to itiTimes (i.e. seconds from start of expt)
    blinkStartTimes = itiTimes + blinkStartTimes/1000;
    blinkEndTimes = itiTimes + blinkEndTimes/1000;

    % check trials don't overlap
    finalBlinkStart=alignRight(blinkStartTimes); % find last saccade in each trial
    if any(finalBlinkStart(1:end-1,end) > itiTimes(2:end))
        disp('final saccade starts later than beginnning of next trial')
        keyboard;
    end

else
    [blinkStartTimes, blinkEndTimes, validBlinks] = deal(NaN); %set to nan, will be removed
end
%% set up [time, code, duration=1];
% saccade start = 52; saccade end = 53;

saccStartCodes = 52 * info.sRT ./ info.sRT; % 52 or NaN
saccEndCodes = 53 * info.sRT ./ info.sRT; % 53 or NaN
blinkStartCodes = 31 * validBlinks; % 31 or NaN
blinkEndCodes = 32 * validBlinks; % 31 or NaN

saccStims = [ col(saccStartTimes), col(saccStartCodes);...
              col(saccEndTimes), col(saccEndCodes);...
              col(blinkStartTimes), col(blinkStartCodes);...
              col(blinkEndTimes), col(blinkEndCodes);...
            ];
saccStims(:,3) = 1; % duration = 1

saccStims = sortrows(saccStims); % order by time
saccStims(isnan(saccStims(:,1)),:) = []; % remove NaN

%% combine with previous and sort again

stimsNew = sortrows( [saccStims; stimsOld] );


%% make "main saccade" = 55
% should take the first saccade after a target appears
% unless they saccaded just before the target (or do something else that
% means there is no saccade detected between target + acceptance) in which
% case it takes the time of acceptance and we will filter these out later

targInds = find(stimsNew(:,2)==48 | stimsNew(:,2)==51); % target onset
acceptInds = find(stimsNew(:,2) == 54); % saccade accepted 
saccInds = find(stimsNew(:,2)==52); % saccade started
if size(targInds,1) ~= size(acceptInds,1)
    warning('different number of target and accepted triggers found');
    keyboard;
end
n = size(targInds,1);
delta = zeros(n,1);
s = zeros(n,3);
for i = 1:length(targInds)
    % find closest saccade to the target
    f = find(saccInds>targInds(i) & saccInds < acceptInds(i),1);
    if isempty(f) % if no saccades between target and accepted, use accepted, will filter out later
        thisSaccInd = acceptInds(i);
    else
        thisSaccInd = saccInds(f);
    end
    % get time difference
    delta(i,1) = stimsNew(thisSaccInd,1) - stimsNew(targInds(i),1);
    s(i,:) = stimsNew(thisSaccInd,:); % get that row
    s(i,2) = 55; % change trigger
    
end

stimsNew = sortrows([stimsNew; s]);

%% check for duplicates

u = unique(stimsNew, 'rows', 'stable');

if size(u,1) ~= size(stimsNew,1)
    disp('some duplicates were detected. deleting them');
    fprintf('%d new triggers added, (%d trials detected)\n', size(u,1) - size(stimsOld,1), size(raw,1));
    stimsNew = u;
end
%%

[a,b] = CountUnique(stimsNew(:,2)); 
disp([b,a])
CheckTrigCodeNumbers(stimsNew(:,2)); % check different triggers match up 

fprintf('is this the correct number of each type of trigger? if not, edit the stims structure. type dbcont to continue\n')
keyboard;
%% examine a few trials to check it all lines up

% targLoc = result.params.targetLocations + result.params.screenSize/2;
% dpp = pix2deg(1, result.params.screenSize(1), 53, 75); % degrees per pixel
% clf;
% for j = 1:25
%     subplot(5,5,j)
%     plot(raw(j,:),'k:')
%     for i=1:sum(~isnan(info.sRT(j,:)))
%         hold on; 
%         plot(raw(j,info.sRT(j,i):(info.sRT(j,i)+info.sDur(j,i))));
%     end
%     DrawCircles(targLoc, 1./dpp)
% end

end