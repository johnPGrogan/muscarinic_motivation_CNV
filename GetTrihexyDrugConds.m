function [isAFirst, ppInds] = GetTrihexyDrugConds(fNames)
% function [conds, ppInds] = GetTrihexyDrugConds()
% Load up the drug conditions from the TrihexyDrugConds.xlsx file, return
% whether each person did A or B on their first session
% 
% Inputs:
%   fNames is a cell array of strings containing TRI_000_0 (ID numbers)

[~,~,raw] = xlsread('../../Experiments/Forms/ACh/TrihexyDrugConds.xlsx');
xlsIDs = raw(:,1);

firstCond = raw(:,3); %get conds

% 1 = A, 0 = B, NaN = neither
isA = NaN(size(firstCond));
isA(strcmp(firstCond,'A')) = 1;
isA(strcmp(firstCond,'B')) = 0;

idStart = regexpi(fNames{1},'Tri_')+4; % first number after TRI_
idEnd = idStart + 4; % first number after TRI_
ids = cellfun(@(x) x(idStart:idEnd), fNames, 'UniformOutput',0);

%% keep only the participants analysed

ppIDs = ids(:,1);
isAFirst = NaN(size(ids,1),1);
for i = 1:size(ppIDs,1)
  ppIDs{i} = sprintf('TRI_%s',ppIDs{i,1}(1:end-2));

  ppInds(i) = max([0 find(strcmp(xlsIDs, ppIDs(i))) ]); % check if IDs match
  
  if ppInds(i)>0 % if so, get the cond
    isAFirst(i,1) = isA(ppInds(i));
  end
end

ppInds(ppInds==0) = NaN;



end