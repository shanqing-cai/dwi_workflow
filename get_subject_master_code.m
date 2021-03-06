function masterCodes = get_subject_master_code(studyIDs, subjIDs, varargin)
%% 
% Optional input arguments: --mat: Use the mat file in lieu of the xls file
%                             Note: scai, use
%                             (SMCG_W510) metaAnalysis/maset_code_xls2mat.m to
%                             convert the xls file into mat.This is mainly
%                             for bypassing the xlsread bug on the BU
%                             servers

%% PYTHON SETTINGS
PYTHON_BIN = '/usr/bin/python';
PY_MASTER_CODE_SCRIPT = 'get_subject_master_code.py';

%% Process options
% bMat = ~isempty(fsic(varargin, '--mat'));

%% Input sanity check
if ~iscell(studyIDs)
    studyIDs = {studyIDs};
end
if ~iscell(subjIDs)
    subjIDs = {subjIDs};
end

if length(studyIDs) ~= length(subjIDs)
    error('Length mismatch between studyIDs and subjIDs');
end

%%
mPath = mfilename('fullpath');
mPath = fileparts(mPath);

check_file(PYTHON_BIN);
pyMCScript = fullfile(mPath, PY_MASTER_CODE_SCRIPT);
check_file(pyMCScript);
cmd = sprintf('%s %s ', PYTHON_BIN, pyMCScript);

for i1 = 1 : numel(studyIDs)
    t_studyID = studyIDs{i1};    
    
    cmd = sprintf('%s%s,', cmd, t_studyID);
end
cmd = sprintf('%s ', cmd(1 : end - 1));
for i1 = 1 : numel(subjIDs)
    t_subjID = subjIDs{i1};    
    
    cmd = sprintf('%s%s,', cmd, t_subjID);
end
cmd = cmd(1 : end - 1);

so = evalc('system(cmd);');

%% Format output
masterCodes = nan(1, length(subjIDs));
res = splitstring(so);
if length(res) ~= length(subjIDs)
   error('Error(s) occurred during python query of SQL database. Command: %s; Output: %s', cmd, so);
end


for i1 = 1 : length(res)
    masterCodes(i1) = str2double(res{i1});
end

return
% 
% %% Load settings from dwi_analysis_settings
% anaSet = load('dwi_analysis_settings.mat');
% 
% %% Locate master code xls file
% assert(isfield(anaSet, 'SUBJECT_MASTER_CODE_FILE'));
% mcfn = anaSet.SUBJECT_MASTER_CODE_FILE;
% 
% if bMat
%     mcfn = strrep(mcfn, '.xls', '.mat');
% end
% 
% check_file(mcfn);
% if ~bMat
%     % Read the xls file    
%     [N, T] = xlsread(mcfn);
% else
%     % Read the mat file
%     load(mcfn);
%     assert(exist('N', 'var') == 1);
%     assert(exist('T', 'var') == 1);
% end
% 
% %% Locate the first line of valid text entry
% row1 = 2;
% while (isempty(T{row1, 1}) || isempty(T{row1, 2}) ...
%        || isempty(T{row1, 3}) || isempty(T{row1, 4}))
%    row1 = row1 + 1;
% end
% 
% 
% %% Iterate through all input studyIDs and subject IDs
% masterCodes = nan(1, length(studyIDs));
% 
% for i1 = 1 : length(studyIDs)
%     t_studyID = studyIDs{i1};
%     t_subjID = subjIDs{i1};
%     
%     bFound = 0;
%     
%     altIDs = get_alt_ids(t_studyID, t_subjID);
%     
%     r = row1;
%     while (r < size(T, 1))
%         xlsID = T{r, 3};
%         
% %         if length(xlsID) > 4 && isequal(xlsID(1 : 4), 'STUT')    % DEBUG
% %             pause(0);
% %         end
% 
%         %--- Build all alternative IDs ---%
% 
%         
%         if ~isempty(fsic(altIDs, xlsID))
%             bFound = 1;
%             break;
%         end
%         
%         r = r + 1;
%     end
%     
%     if bFound
%         mc = strrep(T{r, 2}, 'SL', '');
%     
%         if isnan(str2double(mc))
%             error('Unrecognized master code format in: %s', T{r, 2});
%         end
%     
%         masterCodes(i1) = str2double(mc);
%     else
%         info_log(sprintf('Cannot find the master code for subject %s:%s', ...
%                          t_studyID, t_subjID), ...
%                  '-warn');
%     end
%     
% end
% 
% return
% 
% 
% 
% return

