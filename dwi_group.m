function dwi_group(anaType, grpScheme, varargin)
%% dwi_group
% Input arguments:
%       anaType: analysis type {FA, MD, TNS}
%                TNS - tract density normalized by seed size
%       grpScheme: grouping scheme {byStudy, byGroup}
%
%       --roi roiName: region of interest (ROI) name
%                       (e.g., for FA: lh_vPMC; 
%                              for tns: lh_vMC-lh_vSC)
%       --wm-depth d:   white-matter depth {-1 for gm, 1, 2, 3}
%       --groups grps: groups to analyze
%                      This is compulsory for grpScheme==byGroup
%       --exclude-subjects: exclude specified subjects, seperated by comma
%                           e.g., --exclude-subjects CAT_12,SEQPDS_SEQ02P10
%       --norm-by-all: normalize the FA, MD or other tensor measures by the
%                      average of all ROIs
%       -v | --verbose: verbose mode
%
%%
ANALYSIS_TYPES = {'FA', 'MD', 'TNS'};
GROUPING_SCHEMES = {'byStudy', 'byGroup'};

analysisSettingsMat = 'dwi_analysis_settings.mat';
projInfoMat = 'dwi_project_info.mat';

DEFAULT_PARC = 'aparc12';

%% Visualization options
COLORS = {[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0.5, 1], ...
          [0.5, 1, 0], [0.5, 0, 1], [0.5, 0.5, 0], [0, 1, 1], ...
          [1, 0.25, 0.75], [0.25, 1, 0.75], [0.75, 0.25, 1], [1, 0.75, 0.25], ...
          [0.5, 0.5, 0.25], [0.5, 0.25, 0.75], [0.75, 0.25, 0.5], [0.75, 0.5, 0.5]};
fontSize = 12;

%% Process optional input arguments
args = struct;
args.anaType = anaType;
args.grpScheme = grpScheme;

if ~isempty(fsic(varargin, '--roi'))
    args.roi = varargin{fsic(varargin, '--roi') + 1};
end

if ~isempty(fsic(varargin, '--wm-depth'))
    args.wmDepth = varargin{fsic(varargin, '--wm-depth') + 1};
end

% anaGroups = {};
if ~isempty(fsic(varargin, '--groups'))
    args.groups = varargin{fsic(varargin, '--groups') + 1};
    args.groups = splitstring(args.groups, ',');
end

if ~isempty(fsic(varargin, '--exclude-subjects'))
    args.exclS = varargin{fsic(varargin, '--exclude-subjects') + 1};
    args.exclS = splitstring(args.exclS, ',');
else
    args.exclS = {};
end

bNormByAll = ~isempty(fsic(varargin, '--norm-by-all'));

args.bv = ~isempty(fsic(varargin, '-v')) || ~isempty(fsic(varargin, '--verbose'));

%% Input argument sanity check and processing
if isempty(fsic(ANALYSIS_TYPES, args.anaType))
    error_log(sprintf('Unrecognized anlaysis type: %s', args.anaType));
end
    
if isempty(fsic(GROUPING_SCHEMES, args.grpScheme))
    error_log(sprintf('Unrecognized grouping scheme: %s', args.grpScheme));
end

if isequal(args.anaType, 'TNS')
    assert(length(strfind(args.roi, '-')) == 1);
    rois = splitstring(args.roi, '-');
    seedROI = rois{1};
    targROI = rois{2};

    assert(length(strfind(seedROI, '_')) == 1);
    assert(length(strfind(targROI, '_')) == 1);
    
    items = splitstring(seedROI, '_');
    seedHemi = items{1};
    seedROI = items{2};
    items = splitstring(targROI, '_');
    targHemi = items{1};
    targROI = items{2};
    
    assert(isequal(seedHemi, targHemi)); % May be relaxed in the future for commissural connections
end

%% Check conditionally required input arguments
if isequal(args.anaType, 'FA') || isequal(args.anaType, 'MD')
    if ~isfield(args, 'roi')
        error_log(sprintf('The option --roi for analysis type %s is not supplied', args.anaType));
    end
    
    if ~isfield(args, 'wmDepth')
        error_log(sprintf('The option --wmDepth for analysis type %s is not supplied', args.anaType));
    end
    
end

if isequal(args.grpScheme, 'byGroup') && (~isfield(args, 'group') || length(args.group) == 0)
    error_log(sprintf('--groups not specified under grouping scheme: %s', args.grpScheme));    
end

%% Load dwi_analysis_settings and dwi_project_info 
check_file(analysisSettingsMat);
analysisSettings = load(analysisSettingsMat);

check_file(projInfoMat);
projInfo = load(projInfoMat);

%% Load data from dwi_analysis.py
dat = []; % data
grp = []; % Group info: it can hold things like study ID and group ID

if bNormByAll
    optNormByAll = '--norm-by-all';
else
    optNormByAll = '';
end

for i1 = 1 : size(projInfo.name, 1)
    t_proj = deblank(projInfo.name(i1, :));
    
    ns = length(projInfo.subjIDs{i1}); % Number of subjects
    nsa = 0; % Number of subjects with data available
    
    if args.bv
        info_log(sprintf('Loading data from project %s', t_proj));
    end
    
    for i2 = 1 : size(projInfo.subjIDs{i1}, 1)
        t_sid = deblank(projInfo.subjIDs{i1}(i2, :));
        t_psid = sprintf('%s_%s', t_proj, t_sid);
        
        %-- Exclude subject per input argument --%
        if ~isempty(fsic(args.exclS, t_psid))
            if args.bv
                info_log(sprintf('Excluding subject %s (requested in input arguments)', t_psid));
            end
            continue;
        end
        
        %-- Check group identity of subject --%
        t_grp_id = deblank(projInfo.groupIDs{i1}(i2, :));
        if isfield(args, 'groups') && isempty(fsic(args.groups, t_grp_id))
            if args.bv
                info_log(sprintf('Excluding subject %s (outside groups to be analyzed)', t_psid));
            end
            continue;
        end
        
        %--- Load data ---%
        sDir = fullfile(analysisSettings.DWI_ANALYSIS_DIR, t_psid);
        if isequal(args.anaType, 'FA') || isequal(args.anaType, 'MD')
            tmfile = fullfile(sDir, 'tensor', 'tensor_measures.mat');
            
            if ~isfile(tmfile)
                if args.bv
                    info_log(sprintf('Skipping subject %s, due to missing mat file: %s', ...
                                     t_psid, tmfile), '--warn');
                end
                continue;
            end
            
            t_dat = get_tensor_measure(tmfile, args.anaType, args.roi, args.wmDepth, ...
                                       optNormByAll);
            
            if isnan(t_dat)
                if args.bv
                    info_log(sprintf('Skipping subject %s, due to NaN value in mat file: %s', ...
                                     t_psid, tmfile), '--warn');
                end
                continue;
            end
                        
        elseif isequal(args.anaType, 'TNS')
            connFile = fullfile(sDir, 'conn', ...
                                sprintf('%s_gm_%s.speech.mat', DEFAULT_PARC, seedHemi));
            
            if ~isfile(connFile)
                if args.bv
                    info_log(sprintf('Skipping subject %s, due to missing connectivity mat file: %s', ...
                                     t_psid, connFile), '--warn');
                end
                continue;
            end
            
            t_dat = get_tract_measure(connFile, args.anaType, ...
                                      seedHemi, seedROI, targHemi, targROI);
        else
            error_log(sprintf('Analaysis type %s has not been implemented yet', ...
                              args.anaType));
        end
        
        dat(end + 1) = t_dat;
        
        %--- Add group identity info ---%
        if isequal(args.grpScheme, 'byStudy')
            grp(end + 1) = i1;
        elseif isequal(args.grpScheme, 'byGroup')
            grp(end + 1) = fsic(args.groups, t_grp_id);
        else
            error_log('Unrecognized grouping scheme: %s', args.grpScheme);
        end
        
        nsa = nsa + 1;
    end
    
    if args.bv
        info_log(sprintf('\tData loaded from %d of %d subjects', nsa, ns));
    end

end

%% Some data formatting

%% Statistical analysis and visualization
if isequal(args.grpScheme, 'byStudy')
    ugrps = unique(grp);

    fhdl = figure('Name', sprintf('%s: %s', args.grpScheme, args.anaType), 'Color', 'w');
    set(gca, 'FontSize', fontSize);
    hold on; box on;
    
    [aov1_p, aov1_table] = anova1(dat, grp, 'off');
    for i1 = 1 : numel(ugrps)
        t_grp = ugrps(i1);
        idx = find(grp == t_grp);
        
        if length(idx) > 0
            plot(idx, dat(idx), 'o-', 'Color', COLORS{i1});
            text(idx(1) + 0.5, max(dat(idx)), ...
                 deblank(projInfo.name(t_grp, :)), ...
                 'FontSize', fontSize, 'Color', COLORS{i1});
        end
    end
    
    xs = get(gca, 'XLim');
    ys = get(gca, 'YLim');
    text(xs(1) + 0.05 * range(xs), ys(1) + 0.05 * range(ys), ...
         sprintf('One-way ANOVA: F(%d,%d)=%f; p=%e', ...
                      aov1_table{2, 3}, aov1_table{3, 3}, ...
                      aov1_table{2, 5}, aov1_p), ...
         'fontSize', fontSize);
%     quickText(fhdl, ...
%               sprintf('One-way ANOVA: F(%d,%d)=%f; p=%e', ...
%                       aov1_table{2, 3}, aov1_table{3, 3}, ...
%                       aov1_table{2, 5}, aov1_p), ...
%               fontSize);
    
    xlabel('Subject #');
   	ylabel(strrep(sprintf('%s: %s', args.anaType, args.roi), '_', '\_'));
    
end

return