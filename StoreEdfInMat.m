function StoreEdfInMat(IDs, behDataFolder, exeFile, skipPractice)
% StoreEdfInMat
% convert edf to asc, readAllEdf, store that in result.s

addpath ../../Experiments/MyFuncs/
if exist('IDs') && ~isempty(IDs)
    if ~iscell(IDs)
        IDs = {IDs};
    end
else
    IDs = {'Tri';}; % find all by default
end

if ~exist('behDataFolder','var') || isempty(behDataFolder)
    behDataFolder = 'D:/Backups/Experiments/Results/Trihexy/BehData/';
    if ~exist(behDataFolder, 'file')
      behDataFolder = 'D:/Backups/Trihexy/BehData/';
    end
    if ~exist(behDataFolder, 'file')
        behDataFolder = 'D:/groganj/Work/Trihexy/BehData/';
    end
end

if ~exist('exeFile','var') || isempty(exeFile)
    exeFile = '..\..\Experiments\MyFuncs\@Edf2Mat\private\edf2asc.exe';
    if ~exist(exeFile, 'file')
      exeFile = 'C:\Users\cogneuro\Documents\John\edf_converter\@Edf2Mat\private\edf2asc.exe' ;
    end
    if ~exist(exeFile, 'file')
        error('edf2asc.exe file not found, please specify full path incl file');
    end
end
    
if ~exist('skipPractice','var') || isempty(skipPractice)
    skipPractice = 1;
end


files = what(behDataFolder);


if ~isempty(IDs)
  
  for i = 1:length(IDs)
  %keep only files with this in
    isOut = regexpi(files.mat,IDs{i});
    isOutput(:,i) = ~cellfun(@isempty,isOut);
  end
  isOutputAll = any(isOutput,2);
  
  if skipPractice
      isPrac = ~cellfun(@isempty, regexpi(files.mat, 'practice'));
      isOutputAll(isPrac) = false;
  end
  files.mat = files.mat(isOutputAll);

end
%%

for i = 1:length(files.mat)
  fprintf('\n%s\n',files.mat{i});
  r = load(fullfile(files.path, files.mat{i}));
  result = r.result;

  if isfield(result, 's') % if it has already been processed, check trial numbers match up
      s = result.s;

      % check orig mat file has same num trials
      r1 = load(fullfile(files.path, [result.file(1:8) '.mat']));
      if ~all(length(result.data)== length(result.s) && length(result.data)==length(r1.result.data) )
          fprintf('result had %d trials, %d saccades, and orig has %d trials\n', length(result.data), length(result.s), length(r1.result.data));
          keyboard;
      end
  end
  
  if ~isempty(result.edfFiles) % if there is an edf file
%     if ~isfield(result, 's') % if the edf is not already stored
      edf2ascMat( fullfile(behDataFolder, result.edfFiles), exeFile, 0 ); % convert edf to asc, don't overwrite
      result.s = readAllEDF(result, behDataFolder); % read asc into result.s

%       if any(size(s) ~= size(result.s))
        figure(1);clf; subplot(1,2,1)
        plot(([result.s.saccadeaccepted_t]-[result.s.startITI_t])./1000, ([result.data.saccadeaccepted]-[result.data.startITI]),'x')
        xlabel('sacc RT (s)'); ylabel('beh RT (s)');
        subplot(1,2,2)
        plot(([result.s.saccadeaccepted_t]-[result.s.startITI_t])./1000 - ([result.data.saccadeaccepted]-[result.data.startITI]))
        ylabel('sacc RT - beh RT (s)');
               
        disp(result)
%         result.deleted = 97:104; % indices of deleted trials

%           disp('result.s has changed')
%           save( fullfile(behDataFolder, files.mat{i}), 'result', '-append'); % save
%       end
%     end

    % if it was resumed, copy practice result into new file
    if length(result.edfFiles) > 1
        clear a;
        for j = 1:length(result.edfFiles)
            a(j) = load(fullfile(behDataFolder, [result.edfFiles{j}(1:8) '.mat']));
        end
        
        for j = 1:length(result.edfFiles)
            if isfield(a(j).result, 'practiceResult') && ~isempty(a(j).result.practiceResult)
                result.practiceResult = a(j).result.practiceResult
                keyboard;
%                 save( fullfile(behDataFolder, files.mat{i}), 'result', '-append'); % save
            end
        end
        
    end

    save( fullfile(behDataFolder, files.mat{i}), 'result', '-append'); % save
  end
  
    
        
end

end