function masterCodes = get_subject_master_code(studyIDs, subjIDs, varargin)
%% 
% Optional input arguments: --mat: Use the mat file in lieu of the xls file
%                             Note: scai, use
%                             (SMCG_W510) metaAnalysis/maset_code_xls2mat.m to
%                             convert the xls file into mat.This is mainly
%                             for bypassing the xlsread bug on the BU
%                             servers

%% Config: project name translation table
pnTable = {{'innerspeech', 'InSp'}};


%% Process options
bMat = ~isempty(fsic(varargin, '--mat'));

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

%% Load settings from dwi_analysis_settings
anaSet = load('dwi_analysis_settings.mat');

%% Locate master code xls file
assert(isfield(anaSet, 'SUBJECT_MASTER_CODE_FILE'));
mcfn = anaSet.SUBJECT_MASTER_CODE_FILE;

if bMat
    mcfn = strrep(mcfn, '.xls', '.mat');
end

check_file(mcfn);
if ~bMat
    % Read the xls file    
    [N, T] = xlsread(mcfn);
else
    % Read the mat file
    load(mcfn);
    assert(exist('N', 'var') == 1);
    assert(exist('T', 'var') == 1);
end

%% Locate the first line of valid text entry
row1 = 2;
while (isempty(T{row1, 1}) || isempty(T{row1, 2}) ...
       || isempty(T{row1, 3}) || isempty(T{row1, 4}))
   row1 = row1 + 1;
end


%% Iterate through all input studyIDs and subject IDs
masterCodes = nan(1, length(studyIDs));

for i1 = 1 : length(studyIDs)
    t_studyID = studyIDs{i1};
    t_subjID = subjIDs{i1};
    
    bFound = 0;
    
    r = row1;
    while (r < size(T, 1))
        xlsID = T{r, 3};
        
%         if length(xlsID) > 4 && isequal(xlsID(1 : 4), 'STUT')    % DEBUG
%             pause(0);
%         end

        %--- Build all alternative IDs ---%
        altIDs = cell(1, 2);
        altIDs{1} = [t_studyID, t_subjID];
        altIDs{2} = [t_studyID, '_', t_subjID];
        
        if length(t_subjID) >= 2 && isequal(t_subjID(1), 'S') && ~isnan(str2double(t_subjID(2)))
            bs = 1;
            altIDs{end + 1} = [t_studyID, t_subjID(2 : end)];
            altIDs{end + 1} = [t_studyID, '_', t_subjID(2 : end)];
        else
            bs = 0;
        end
        
        if length(t_subjID) > 6 && isequal(t_subjID(1 : 3), 'SEQ') && ...
                (isequal(t_subjID(6), 'P') || isequal(t_subjID(6), 'C'))
            bSEQPDS = 1;
            altIDs{end + 1} = strrep(t_subjID, 'SEQ', 'SEQPDS');            
        else
            bSEQPDS = 0;
        end
        
        if length(t_subjID) > 4 && isequal(t_subjID(1 : 3), 'FRS') && ...
                ~isnan(str2double(t_subjID(4 : end)))
            bFRS = 1;
            altIDs{end + 1} = t_subjID;
        else
            bFRS = 0;
        end
        
        if length(t_subjID) > 5 && isequal(t_subjID(1 : 4), 'CCRS') && ...
                ~isnan(str2double(t_subjID(5 : end)))
            bCCRS = 1;
            altIDs{end + 1} = t_subjID;
        else
            bCCRS = 0;
        end

        if length(t_subjID) > 4 && isequal(t_subjID(1 : 3), 'SEQ') && ...
                ~isnan(str2double(t_subjID(4 : end)))
            bSEQ = 1;
            altIDs{end + 1} = t_subjID;
        else
            bSEQ = 0;
        end
        
        %--- Find alternative project names from the pnTable ---% 
        altStudyIDs = get_alt_study_ids(pnTable, t_studyID);
        if ~isempty(altStudyIDs)
            for i2 = 1 : numel(altStudyIDs)
                altIDs{end + 1} = [altStudyIDs{i2}, t_subjID];
                altIDs{end + 1} = [altStudyIDs{i2}, '_', t_subjID];
                if bs
                    altIDs{end + 1} = [altStudyIDs{i2}, t_subjID(2 : end)];
                    altIDs{end + 1} = [altStudyIDs{i2}, '_', t_subjID(2 : end)];
                end
            end
        end
        
        if ~isempty(fsic(altIDs, xlsID))
            bFound = 1;
            break;
        end
        
        r = r + 1;
    end
    
    if bFound
        mc = strrep(T{r, 2}, 'SL', '');
    
        if isnan(str2double(mc))
            error('Unrecognized master code format in: %s', T{r, 2});
        end
    
        masterCodes(i1) = str2double(mc);
    else
        info_log(sprintf('Cannot find the master code for subject %s:%s', ...
                         t_studyID, t_subjID), ...
                 '-warn');
    end
    
end

return

%% Subroutines
function altStudyIDs = get_alt_study_ids(pnTable, t_studyID)
altStudyIDs = {};
for i1 = 1 : length(pnTable)
    if isequal(pnTable{i1}{1}, t_studyID)
        altStudyIDs = pnTable{i1}(2 : end);
    end
end

return

