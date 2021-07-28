function stims = InsertTriggerCodes(trigChannel, sampleTime, trigDeltaThresh, ET, triggerScaling)
% stims = InsertTriggerCodes(trigChannel, sampleTime, trigDeltaThresh)
% finds trigger values based on threshold for diff()
% gets the delay between on&off as the trigger value, rounds to nearest
% integer
% Inputs:   trigChannel is a vector of voltages recorded from the channel
%               with trigger codes as a voltage step where the duration (ms)
%               corresponds to the value of the trigger code
%           sampleTime is a vector of times (ms) of each sample taken
%           trigDeltaThresh is the threshold for successive same triggers
%               being removed (default = 1, RPM should be < 0.3).
% 
% returns stims = [time, code, duration] of each trigger 
% 
% John Grogan, 2019.
if ~exist('triggerScaling', 'var'); triggerScaling = 3; end

%% get times of triggers

trigInds = findregions(trigChannel==254);

trigTimes = sampleTime(trigInds); % convert into times of onset + offset

trigLength = diff(trigTimes,1,2) * 1000; % length of code in ms

trigCodes = round(trigLength/triggerScaling)*triggerScaling; % round to get code values
% codes are 3:3:36, there is some variance in timing due to device drift,
% so we can correct by doing round(trigLength/3)*3


%% use ET.events as trigCodes - avoids rounding errors
% we are only taking the trigValues, not the timings, the times are still
% from the onset of the triggers

if exist('ET','var') && ~isempty(ET)
    if length(trigCodes) ~= length(ET.event)
      disp('trigCodes and ET.event do not have the same number of trigger events')
      keyboard;
    else
      t = [trigCodes, ET.event(:,2)];
      disp(t(trigCodes~=ET.event(:,2),:))
      a = input('accept ET.event values as trigger values? 1/0: ');
      if a==1
        trigCodes = ET.event(:,2);
      else
        disp('you have chosen not to use ET.event as the trigger codes. type dbcont to continue')
        keyboard;
      end
    end
end
%% insert trig codes back into data at correct time stamp

stims(:,1) = trigTimes(:,1); % time in ms
stims(:,2) = trigCodes; % ID of code
stims(:,3) = trigLength; % duration in ms

%% some codes were sent before&after an event, so the event lies between the end of one and start of another 

sameTrig = diff(trigCodes)==0; % get when code was repeated
trigDelta = trigTimes(2:end,1) - trigTimes(1:end-1,2); % time between one ending and another starting

if ~exist('trigDeltaThresh','var')
    trigDeltaThresh = 1;
end
removeDupl = trigDelta < trigDeltaThresh & sameTrig;

if any(removeDupl)
  warning('removing some duplicated trigger events');
  stims(removeDupl,:) = [];
end

%% count how many trigger codes there are - adjust if necessary
% count how many trigger codes there are 
[a,b] = CountUnique(stims(:,2));
disp([b, a])

CheckTrigCodeNumbers(stims(:,2)); % check things add up
% f = find(stims(2:end,2)==21 & stims(1:end-1,2)>=21)
% stims(f+1,:) = [];

fprintf('is this the correct number of each type of trigger? if not, edit the stims structure. type dbcont to continue\n')
keyboard;
end