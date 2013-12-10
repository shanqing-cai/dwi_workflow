function ys = get_subject_demo_info(infoType, varargin)
%%
% infoType: {age}

%% CONFIG
demoXlsFN = '/speechlab/2/jtour/SID/SL_demos.xls';
demoMatFN = strrep(demoXlsFN, '.xls', '.mat');

% WANRING: Ad hoc! May likely be screwed up by changes in the xls file
xlsAgeCol = 8;
xlsAgeRowOffset = -2;
xlsMasterCodeCol = 1;
xlsSIDCol = 2;

%%
if isnumeric(varargin{1}(1))
    inputMode = 'masterCode';
    masterCodes = varargin{1};
    
    ys = nan(size(masterCodes));
else
    inputMode = 'studySubjIDs';
    studyIDs = varargin{1};
    subjIDs = varargin{2};
    
    assert(iscell(studyIDs));
    assert(iscell(subjIDs));
    
    assert(length(studyIDs) == length(subjIDs));
    ys = nan(size(studyIDs));
end

%%
bXls = ~isempty(fsic(varargin, '--xls'));
bMat = ~isempty(fsic(varargin, '--mat'));

if bXls && bMat
    error('Incompatible options --xls and --mat are used simultaneously');
end

%% Load info from xls file
if bXls
    check_file(demoXlsFN);
    [N, T] = xlsread(demoXlsFN);
elseif bMat
    check_file(demoMatFN);
    load(demoMatFN);
    
    assert(exist('N', 'var') == 1);
    assert(exist('T', 'var') == 1);
end

%%

if ~bXls && ~bMat
    tpath = which('db_search');
    if isempty(tpath)
        error('Cannot find path to function db_search');
    end

    if isequal(infoType, 'age')
        r = db_search('subject.id', 'subject.age');
    end

    for i1 = 1 : length(masterCodes)
        idx = find(r{1} == masterCodes(i1));
        if length(idx) == 0
            info_log(sprintf('Cannot find the %s info for subject masterCode == %d', ...
                             infoType, masterCodes(i1)), '-warn');
            continue;
        end
        if length(idx) > 1
            info_log(sprintf('More than one entries found for subject masterCode == %d. Will use the first entry', ...
                             infoType, masterCodes(i1)), '-warn');
            idx = idx(1);    
        end

        ys(i1) = r{2}(idx);
    end
else
    % Search for the subject
    if isequal(inputMode, 'masterCode')
        masterCodeCol = T(:, xlsMasterCodeCol);
        
        for i1 = 1 : length(masterCodes)
            t_mc = masterCodes(i1);
            t_mcStr = sprintf('SL%.4d', t_mc);

            idx = fsic(masterCodeCol, t_mcStr);
            if length(idx) == 0
                info_log(sprintf('Cannot find the %s info for subject masterCode == %d', ...
                                 infoType, masterCodes(i1)), '-warn');
                continue;
            end

            if length(idx) > 1
                info_log(sprintf('More than one entries found for subject masterCode == %d. Will use the first entry', ...
                                 infoType, masterCodes(i1)), '-warn');
                idx = idx(1);
            end

            ys(i1) = N(idx + xlsAgeRowOffset, xlsAgeCol);
        end
    elseif isequal(inputMode, 'studySubjIDs')
        for i1 = 1 : length(studyIDs)
            sidCol = T(:, xlsSIDCol);
            
            t_studyID = studyIDs{i1};
            t_subjID = subjIDs{i1};
            
            altIDs = get_alt_ids(t_studyID, t_subjID);

            idx = [];
            for i2 = 1 : numel(altIDs)
                idx = [idx, fsic(sidCol, altIDs{i2})];
            end
            
            if length(idx) == 0
                info_log(sprintf('Cannot find the %s info for subject proj == %s; sid = %s', ...
                                 infoType, t_studyID, t_subjID), '-warn');
                continue;
            end

            if length(idx) > 1
                info_log(sprintf('More than one entries found for subject proj == %s; sid = %s. Will use the first entry', ...
                                 t_studyID, t_subjID), '-warn');
                idx = idx(1);
            end

            ys(i1) = N(idx + xlsAgeRowOffset, xlsAgeCol);
        end
    end
        
end
return