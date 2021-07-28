function [EEG, epochFlags] = GetEpochFlags(EEG, updateManual)
% [EEG, epochFlags] = GetEpochFlags(EEG, updateManual)
% converts the binary flagged events within each epoch into indices of
% flags for each trial (across all events within an epoch). This is stored
% in EEG.epoch.trailFlags. 
% If updateManual = 1, then it will force those to overwrite the
% epochFlags, such that any unmarked will be set all to zero, and any
% marked will have the 1st flag only set to 1. Default is 0.
% 

%% convert flag codes from binary into indices 1:16
ntrial  = EEG.trials;
nbin    = EEG.EVENTLIST.nbin;
oldflag = zeros(1,ntrial);

for b=1:nbin
    for i=1:ntrial
        if length(EEG.epoch(i).eventlatency) == 1
            binix = [EEG.epoch(i).eventbini];
            if iscell(binix)
                binix = cell2mat(binix);
            end
            if ismember(b, binix)
                flagx = [EEG.epoch(i).eventflag];
                if iscell(flagx)
                    flagx = cell2mat(flagx);
                end
                oldflag(b,i)   = flagx;
            else
                oldflag(b,i) =0;
            end
        elseif length(EEG.epoch(i).eventlatency) > 1
            indxtimelock = find(cell2mat(EEG.epoch(i).eventlatency) == 0,1,'first'); % catch zero-time locked type
            if ismember(b, EEG.epoch(i).eventbini{indxtimelock})
                oldflag(b,i)   = EEG.epoch(i).eventflag{indxtimelock};
            else
                oldflag(b,i) =0;
            end
        else
            errorm  = 1;
        end
    end
end

%% get which flags they were
flagbit    = bitshift(1, 0:15);

for b=1:nbin
    for j=1:16
        C(:,j,b) = bitand(flagbit(j), oldflag(b,:));
    end
end

epochFlags = sum(C>0, 3); % sum across bins

%% overwrite by manual rejection (from eegplot)?

if exist('updateManual', 'var') && updateManual
    rejManual = EEG.reject.rejmanual;
    if ~isempty(rejManual)
        epochFlags(rejManual==0,:) = 0;
        epochFlags(rejManual==1,1) = 1;
    end
end


%% store

for i = 1:ntrial
    EEG.epoch(i).trialflag = epochFlags(i,:);
end

    