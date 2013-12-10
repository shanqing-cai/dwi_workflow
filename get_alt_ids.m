function altIDs = get_alt_ids(t_studyID, t_subjID)
%% Config: project name translation table
pnTable = {{'innerspeech', 'InSp'}};

SL_DEMO_MAT_FN = '/speechlab/2/jtour/SID/SL_demos.mat';
USE_MAT_STUDY_CODE_PROJS = {'SDAP'};
STUDY_CODE_COL = 2;
DATA_DIR_NAME_COL = 3;

%%
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

if ~isempty(fsic(USE_MAT_STUDY_CODE_PROJS, t_studyID))
    load(SL_DEMO_MAT_FN);
    assert(exist('N', 'var') == 1);
    assert(exist('T', 'var') == 1);
    
    bFound = 0;
    for i1 = 1 : size(T, 1)
        if ~isempty(fsic(altIDs, T{i1, DATA_DIR_NAME_COL}))
            bFound = 1;
            break;
        end
    end
    
    if bFound
        altIDs{end + 1} = T{i1, STUDY_CODE_COL};
        altIDs{end + 1} = [t_studyID, '_', T{i1, STUDY_CODE_COL}];
    end
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
        
return

%% Subroutines
function altStudyIDs = get_alt_study_ids(pnTable, t_studyID)
altStudyIDs = {};
for i1 = 1 : length(pnTable)
    if isequal(pnTable{i1}{1}, t_studyID)
        altStudyIDs = pnTable{i1}(2 : end);
    end
end