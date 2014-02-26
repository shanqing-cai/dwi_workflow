function varargout = dwi_group(anaType, grpScheme, varargin)
%% dwi_group
% Input arguments:
%       anaType: analysis type {FA, MD, TNS}
%                TNS - tract density normalized by seed size
%       grpScheme: grouping scheme {byStudy, byGroup, rep, byAge}
%                   rep - Reproducibility test: multiple sessions on the
%                         same subject
%                   age - Analayze age effects
%
%       --projs: limit the analysis to specific projects: string list of
%                projects, separeted by commas (e.g., 'STUT,RHY')
%       --roi roiName: region of interest (ROI) name
%                       (e.g., for FA: lh_vPMC; 
%                              for tns: lh_vMC-lh_vSC)
%       --wm-depth d:   white-matter depth {-1 for gm, 1, 2, 3}
%       --path-name pathName: path name to be analyzed under the "path" 
%                             analysis mode, e.g., lh_SLFT
%       --path-meas pathMeasure: measure of path to be analyzed under the
%                                "path" mode, e.g., FA_Avg_Weight
%       --groups grps: groups to analyze
%                      This is compulsory for grpScheme==byGroup
%       --exclude-subjects: exclude specified subjects, seperated by comma
%                           e.g., --exclude-subjects CAT_12,SEQPDS_SEQ02P10
%       --norm-by-all: normalize the FA, MD or other tensor measures by the
%                      average of all ROIs
%       --no-self:     remove the self-projections in the TNS analysis
%                      Self-projections refers to the projection from an ROI to itself                    
%       --rand-other: Use random other individual subject to compare with
%                     the TNS measure of the repeated subject.
%       --symm-tract: Use symmetric tracts for tractography results 
%                     (i.e., average the two symmetrical elements of the
%                     connectivity matrix)
%       -v | --verbose: verbose mode
%
%%
ANALYSIS_TYPES = {'FA', 'MD', 'TNS', 'path'};
GROUPING_SCHEMES = {'byStudy', 'byGroup', 'rep', 'byAge'};

analysisSettingsMat = 'dwi_analysis_settings.mat';
projInfoMat = 'dwi_project_info.mat';

DEFAULT_PARC = 'aparc12';

%% Other constants
MIN_TNS = 1e-6;

%--- Visualization settings ---%
fontSize = 12;

FA_LIMS = [0.1, 0.6];

%% Visualization options
COLORS = {[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0.5, 1], ...
          [0.5, 0.25, 0], [0.5, 0, 1], [0.5, 0.5, 0], [0, 1, 1], ...
          [1, 0.25, 0.75], [0.25, 1, 0.75], [0.75, 0.25, 1], [1, 0.75, 0.25], ...
          [0.5, 0.5, 0.25], [0.5, 0.25, 0.75], [0.75, 0.25, 0.5], [0.75, 0.5, 0.5]};

%-- Colors for group comparisons --%
G_COLORS.nrm = [0, 0, 1];
G_COLORS.pds = [1, 0, 0]; % PDS: persistent developmental stuttering
G_COLORS.sdp = [1, 0.5, 0]; % SDP: spasmodic dysphonia

%% Process optional input arguments
args = struct;
args.anaType = anaType;
args.grpScheme = grpScheme;

if ~isempty(fsic(varargin, '--proj'))
    args.projs = varargin{fsic(varargin, '--proj') + 1};
    args.projs = splitstring(args.projs, ',');
end

if ~isempty(fsic(varargin, '--roi'))
    args.roi = varargin{fsic(varargin, '--roi') + 1};
end

if ~isempty(fsic(varargin, '--wm-depth'))
    args.wmDepth = varargin{fsic(varargin, '--wm-depth') + 1};
end

if ~isempty(fsic(varargin, '--path-name'))
    args.pathName = varargin{fsic(varargin, '--path-name') + 1};
end

if ~isempty(fsic(varargin, '--path-meas'))
    args.pathMeas = varargin{fsic(varargin, '--path-meas') + 1};
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

args.bNormByAll = ~isempty(fsic(varargin, '--norm-by-all'));
args.bNoSelf = ~isempty(fsic(varargin, '--no-self'));

args.bRandOther = ~isempty(fsic(varargin, '--rand-other'));

args.symmTract = ~isempty(fsic(varargin, '--symm-tract'));

args.bv = ~isempty(fsic(varargin, '-v')) || ~isempty(fsic(varargin, '--verbose'));

%% Input argument sanity check and processing
if isempty(fsic(ANALYSIS_TYPES, args.anaType))
    error_log(sprintf('Unrecognized anlaysis type: %s', args.anaType));
end
    
if isempty(fsic(GROUPING_SCHEMES, args.grpScheme))
    error_log(sprintf('Unrecognized grouping scheme: %s', args.grpScheme));
end

if isequal(args.anaType, 'TNS')
    if ~(isequal(args.roi, 'lh_all') || isequal(args.roi, 'rh_all'))
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

end


%% Check conditionally required input arguments
if isequal(args.anaType, 'FA') || isequal(args.anaType, 'MD')
    if ~isfield(args, 'roi')
        error_log(sprintf('The option --roi for analysis type %s is not supplied', args.anaType));
    end
    
    if ~isfield(args, 'wmDepth')
        error_log(sprintf('The option --wm-depth for analysis type %s is not supplied', args.anaType));
    end
elseif isequal(args.anaType, 'path')
    if ~isfield(args, 'pathName')
        error_log('Under analysis type "path", the compulsory input argument --path-name is missing');
    end
    if ~isfield(args, 'pathMeas')
        error_log('Under analysis type "path", the compulsory input argument --path-meas is missing');
    end
end

if isequal(args.grpScheme, 'byGroup') && (~isfield(args, 'groups') || isempty(args.groups))
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
projs = {};     % Project names
sids = {};      % Subject IDs


if args.bNormByAll
    optNormByAll = '--norm-by-all';
else
    optNormByAll = '';
end

prev_roiNames = {};

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
            
            if ~isequal(args.roi, 'all')
                t_dat = get_tensor_measure(tmfile, args.anaType, args.roi, args.wmDepth, ...
                                           optNormByAll);
            else
                [t_dat, t_roiNames] = get_tensor_measure(tmfile, args.anaType, args.roi, args.wmDepth, ...
                                                         optNormByAll);
                if ~isempty(prev_roiNames)
                    assert(isequal(t_roiNames, prev_roiNames));
                end
            end
            
            
            if ~isequal(args.roi, 'all') 
                if isnan(t_dat)
                    if args.bv
                        info_log(sprintf('Skipping subject %s, due to NaN value in mat file: %s', ...
                                         t_psid, tmfile), '--warn');
                    end
                    
                    continue;
                end                
            end
                        
        elseif isequal(args.anaType, 'TNS')
            if isequal(args.roi, 'lh_all') || isequal(args.roi, 'rh_all')
                seedHemi = args.roi(1 : 2);
            end
            connFile = fullfile(sDir, 'conn', ...
                                sprintf('%s_gm_%s.speech.mat', DEFAULT_PARC, seedHemi));
            
            if ~isfile(connFile)
                if args.bv
                    info_log(sprintf('Skipping subject %s, due to missing connectivity mat file: %s', ...
                                     t_psid, connFile), '--warn');
                end
                continue;
            end
            
            if ~(isequal(args.roi, 'lh_all') ||     isequal(args.roi, 'rh_all'))
                t_dat = get_tract_measure(connFile, args.anaType, ...
                                          seedHemi, seedROI, targHemi, targROI);
                                      
                if args.symmTract
                    t_dat_sym = get_tract_measure(connFile, args.anaType, ...
                                                  targHemi, targROI, seedHemi, seedROI);
                    t_dat = (t_dat + t_dat_sym) / 2;
                end
                
            else
                [t_dat, t_roiNames] = get_tract_measure(connFile, args.anaType, ...
                                                        seedHemi, 'all');
                if ~isempty(prev_roiNames)
                    assert(isequal(t_roiNames, prev_roiNames));
                end
            end
        elseif isequal(args.anaType, 'path')
            dpathDataMat = fullfile(sDir, 'dpath', 'dpath_data.mat');
            if ~isfile(dpathDataMat)
                if args.bv
                    
                    info_log(sprintf('Skipping subject %s, due to missing dpath data mat file: %s', ...
                                     t_psid, dpathDataMat), '--warn');
                end
                continue;
            end
            
            pathData = load(dpathDataMat);
            
            %-- Locate the path --%
            allPaths = fields(pathData);
            idxMatch = strmatch(lower(args.pathName), allPaths);
            if length(idxMatch) ~= 1
                error_log(sprintf('Cannot find exactly one entry in path data matching the path name: %s', args.pathName));
            end
            fld = allPaths{idxMatch};
            
            t_dat = pathData.(fld).(args.pathMeas);
            
        else
            error_log(sprintf('Analaysis type %s has not been implemented yet', ...
                              args.anaType));
        end
        
        if isfield(args, 'roi') && ...
           (isequal(args.roi, 'all') || isequal(args.roi, 'lh_all') || isequal(args.roi, 'rh_all'))
            if isvector(t_dat)
                dat = [dat, t_dat];
            else
                dat = cat(3, dat, t_dat);
            end
            prev_roiNames = t_roiNames;
        else
            dat(end + 1) = t_dat;
        end
        
        %--- Add project and subject ID info ---%
        projs{end + 1} = t_proj;
        sids{end + 1} = t_sid;
        
        %--- Add group identity info ---%
        if isequal(args.grpScheme, 'byStudy')
            grp(end + 1) = i1;
        elseif isequal(args.grpScheme, 'byGroup')
            grp(end + 1) = fsic(args.groups, t_grp_id);
        elseif isequal(args.grpScheme, 'rep')            
            grp(end + 1) = NaN;
        elseif isequal(args.grpScheme, 'byAge')
            grp(end + 1) = i1;
        else
            error_log('Unrecognized grouping scheme: %s', args.grpScheme);
        end
        
        nsa = nsa + 1;
    end
    
    if args.bv
        info_log(sprintf('\tData loaded from %d of %d subjects', nsa, ns));
    end

end

%% Get the master code of all subjects
if args.bv
    info_log(sprintf('Getting the master code of all (%d) subjects...', length(projs)));
end

masterCodes = get_subject_master_code(projs, sids, '--mat');
if args.bv
    info_log(sprintf('Done.'));
end

%% Get a list of all the repeats (only two in each row)
if isequal(args.grpScheme, 'rep')
    repTab = nan(0, 2);
    repProjs = cell(0, 2);
    u_masterCodes = unique(masterCodes(~isnan(masterCodes)));
    for i1 = 1 : length(u_masterCodes)
        idxs = find(masterCodes == u_masterCodes(i1));
        if length(idxs) > 1
            for i2 = 2 : length(idxs)
                repTab = [repTab; [idxs(i2 - 1), idxs(i2)]];
                repProjs = [repProjs; {projs{idxs(i2 - 1)}, projs{idxs(i2)}}];
            end
        end    
    end
    if args.bv
        info_log(sprintf('Rep analysis: found %d pairs of repeated experiments on same subject', ...
                         size(repTab, 1)));
    end
end

%% Age statistics
if isequal(args.grpScheme, 'byAge')
    ages = get_subject_demo_info('age', projs, sids, '--mat');
end

%% Select projects
if isfield(args, 'projs') && ~isempty(args.projs)
    idxPres = [];
    
    for i1 = 1 : length(args.projs)
        idxPres = [idxPres, fsic(projs, args.projs{i1})];
    end
    
    if isempty(idxPres)
        error_log('There are no subjects from the projects: %s', ...
                  varargin{fsic(varargin, '--projs') + 1});
    end
    
    projs = projs(idxPres);
    dat = dat(idxPres);
    grp = grp(idxPres);
    sids = sids(idxPres);
    
    if exist('masterCodes', 'var')
        masterCodes = masterCodes(idxPres);
    end
    if exist('ages', 'var')
        ages = ages(idxPres);
    end
end

%% Translate projs to projids (numerical)
uprojs = unique(projs);
projids = nan(size(projs));
for i1 = 1 : numel(projs)
    projids(i1) = fsic(uprojs, projs{i1});
end

%% Statistical analysis and visualization
if isfield(args, 'roi')
    measName = strrep(sprintf('%s: %s', args.anaType, args.roi), '_', '\_');
elseif isfield(args, 'pathName')
    measName = strrep(sprintf('%s: %s - %s', args.anaType, args.pathName, args.pathMeas), '_', '\_');
end
    
if args.bNormByAll
    measName = [measName, ' (norm. by allAvg)'];
end

if isequal(args.grpScheme, 'byStudy') || isequal(args.grpScheme, 'byAge')
    ugrps = unique(grp);

    if isequal(args.grpScheme, 'byAge')
        fhdl = figure('Name', sprintf('%s: %s', args.grpScheme, args.anaType), 'Color', 'w', ...
                      'Position', [100, 200, 1200, 500]);
    else
        fhdl = figure('Name', sprintf('%s: %s', args.grpScheme, args.anaType), 'Color', 'w');  
        set(gca, 'FontSize', fontSize);
        hold on; box on;
    end

    
    if isequal(args.grpScheme, 'byAge')
        spN = 2;
        spM = ceil(numel(ugrps) / 2);
    end
    
    if length(size(dat)) == 2 && isvector(dat) %--- Scalar measure from each subject ---%
        [aov1_p, aov1_table] = anova1(dat, grp, 'off');
        
        for i1 = 1 : numel(ugrps)
            if isequal(args.grpScheme, 'byAge')
                subplot(spN, spM, i1);
                set(gca, 'FontSize', fontSize);
                hold on; box on;
            end
            
            t_grp = ugrps(i1);
            idx = find(grp == t_grp);
            
            if isequal(args.grpScheme, 'byStudy')
                plot(idx, dat(idx), 'o-', 'Color', COLORS{i1});
            else
                plot(ages(idx), dat(idx), 'o', 'Color', COLORS{i1});
            end

            if isequal(args.grpScheme, 'byAge')
                if  isequal(args.anaType, 'FA')
                    set(gca ,'YLim', FA_LIMS);
                end
                title(sprintf('%s: %s', deblank(projInfo.name(t_grp, :)), args.anaType));
                
                %--- Peform linear correlation ---%
                [k_lin, r2_lin, p_lin] = lincorr(ages(idx), dat(idx));
                xs = get(gca, 'XLim');
                ys = get(gca, 'YLim');
                plot(xs, k_lin(1) + k_lin(2) * xs, '--', 'Color', COLORS{i1});

                r_lin = ((k_lin(2) > 0) * 2 - 1) * sqrt(r2_lin);
                
                text(xs(1) + 0.05 * range(xs), ys(2) - 0.06 * range(ys), ...
                     sprintf('Lincorr: R=%.3f, p=%.3f', r_lin, p_lin), ...
                     'FontSize', fontSize);
                 
                %--- Spearman's correlation ---%
                [rho_sp, t_sp, p_sp] = spear(ages(idx).', dat(idx).');
                text(xs(1) + 0.05 * range(xs), ys(2) - 0.14 * range(ys), ...
                     sprintf('Spearman: rho=%.3f, p=%.3f', rho_sp, p_sp), ...
                     'FontSize', fontSize);
            else
                if length(idx) > 0
                    if isequal(args.grpScheme, 'byStudy')
                        text(idx(1) + 0.5, max(dat(idx)), ...
                             deblank(projInfo.name(t_grp, :)), ...
                             'FontSize', fontSize, 'Color', COLORS{i1});
                    elseif isequal(args.grpScheme, 'byAge')
                        text(ages(idx(1)), dat(idx(1)), ...
                             deblank(projInfo.name(t_grp, :)), ...
                             'FontSize', fontSize, 'Color', COLORS{i1});
                    end
                end
            end
        end

        xs = get(gca, 'XLim');
        ys = get(gca, 'YLim');
        text(xs(1) + 0.05 * range(xs), ys(1) + 0.05 * range(ys), ...
             sprintf('One-way ANOVA: F(%d,%d)=%f; p=%e', ...
                          aov1_table{2, 3}, aov1_table{3, 3}, ...
                          aov1_table{2, 5}, aov1_p), ...
             'fontSize', fontSize);

        if isequal(args.grpScheme, 'byStudy')
            xlabel('Subject #');
        else
            xlabel('Age (y.o.)');
        end

        
    elseif length(size(dat)) == 2 %-- Vector measure from each subject --%
        t_legend = {};
        for i1 = 1 : numel(ugrps)
            t_grp = ugrps(i1);
            idx = find(grp == t_grp);

            if length(idx) > 0
                mean_dat = nanmean(dat(:, idx), 2);
                if isequal(args.grpScheme, 'byStudy')
                    plot(1 : length(prev_roiNames), mean_dat, '-', 'Color', COLORS{i1});
                    t_legend{end + 1} = deblank(projInfo.name(t_grp, :));
                elseif isequal(args.grpScheme, 'byAge')
                    error('grpScheme=byAge currently not supported under roi=all');
                end
            end
        end
        
        xlabel('ROI');        
    elseif length(size(dat)) == 3 %-- Matrix measure from each subject --%        
        t_legend = {};
        for i1 = 1 : numel(ugrps)
            t_grp = ugrps(i1);
            idx = find(grp == t_grp);

            if length(idx) > 0
                mean_mat = nanmean(dat(:, :, idx), 3);
                
                if isequal(args.anaType, 'TNS') && args.bNoSelf
                    for k1 = 1 : size(mean_mat, 1)
                        mean_mat(k1, k1) = NaN;
                    end
                end
                
                mean_dat = reshape(mean_mat, size(dat, 1) * size(dat, 2), 1);
                mean_dat = mean_dat(~isnan(mean_dat));
                
                if isequal(args.grpScheme, 'byStudy')
                    plot(1 : length(mean_dat), mean_dat, '-', 'Color', COLORS{i1});
                    t_legend{end + 1} = deblank(projInfo.name(t_grp, :));
                elseif isequal(args.grpScheme, 'byAge')
                    error('grpScheme=byAge currently not supported under roi=all');
                end
            end
        end
        
        set(gca, 'YScale', 'log');
        
        fpos = get(fhdl, 'Position');
        fpos(3) = fpos(3) * 2;
        set(fhdl, 'Position', fpos);
        
        xlabel('ROI-ROI pair');
    end
    
    ylabel(measName);
    
    if length(size(dat)) >= 2 && ~isvector(dat)
        legend(t_legend);        
    end
   
elseif isequal(args.grpScheme, 'rep')
    if (isequal(args.anaType, 'FA') || isequal(args.anaType, 'MD')) && ~isvector(dat)
        figure('Position', [100, 200, 1400, 500], 'Color', 'w');
        spN = 2;
        spM = ceil(size(repTab, 1) / 2);

        rhos_self = nan(size(repTab, 1), 1);
        rhos_so = nan(size(repTab, 1), 1);
        r2s_self = nan(size(repTab, 1), 1);
        r2s_so = nan(size(repTab, 1), 1);
        for i1 = 1 : size(repTab, 1)
            idxOther = setxor(1 : size(dat, 2), repTab(i1, :));

            dat_rx = dat(:, repTab(i1, 1));
            dat_ry = dat(:, repTab(i1, 2));
            dat_ro = nanmean(dat(:, idxOther), 2);

            subplot(spN, spM, i1);
            hold on;
            plot(dat_rx, dat_ry, 'bo');
            plot(dat_rx, dat_ro, 'ro');

            xs = get(gca, 'XLim');
            ys = get(gca, 'YLim');
            lims = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
            plot(lims, lims, '-', 'Color', [0.5, 0.5, 0.5]);
            set(gca, 'XLim', lims, 'YLim', lims);
            grid on;

            [rho_self, t_self, p_sp_self] = spear(dat_rx, dat_ry);
            [rho_so, t_so, p_sp_so] = spear(dat_rx, dat_ro);

            [k_self, r2_self, p_lin_self] = lincorr(dat_rx, dat_ry);
            [k_so, r2_so, p_lin_so] = lincorr(dat_rx, dat_ro);

            rhos_self(i1) = rho_self;
            rhos_so(i1) = rho_so;
            r2s_self(i1) = r2_self;
            r2s_so(i1) = r2_so;

            text(xs(1) + 0.01 * range(xs), ys(2) - 0.06 * range(ys), ...
                 sprintf('Self-self: rho=%.3f, R^2=%.3f', rho_self, r2_self), ...
                 'Color', 'b');
            text(xs(1) + 0.01 * range(xs), ys(2) - 0.12 * range(ys), ...
                 sprintf('Self-other: rho=%.3f, R^2=%.3f', rho_so, r2_so), ...
                 'Color', 'r');

            xlabel(sprintf('measName from Session 1 (%s)', projs{repTab(i1, 1)}));
            ylabel(sprintf('measName from Session 2 (%s)', projs{repTab(i1, 2)}));
        end

        figure('Color', 'w'); 
        set(gca, 'FontSize', fontSize);
        hold on;
        plot(r2s_self, r2s_so, 'ko');
        lims = [0, 1];
        set(gca, 'XLim', lims, 'YLim', lims);
        grid on;
        plot(lims, lims, '-', 'Color', [0.5, 0.5, 0.5]);
        xlabel('R^2 of self-self correlation');
        ylabel('R^2 of self-other correlation');  
%         end        
    else
        bPlot = (nargout == 0);
        
        if length(size(dat)) == 2   % Single connection
            dat_rx = dat(repTab(:, 1));
            dat_ry = dat(repTab(:, 2));
            
            %-- Perform linear correlation --%
            [lc_k, lc_r2, lc_p] = lincorr(dat_rx, dat_ry);
            lc_r = sqrt(lc_r2) * sign(lc_k(2));
            
            %-- Perform Spearman's correlation --%
            [sp_r, sp_t, sp_p] = spear(dat_rx(:), dat_ry(:));

            if ~bPlot
                varargout{1} = [lc_r, sp_r];
                return
            else
                figure;
                hold on;
                for i1 = 1 : length(dat_rx) 
                        plot(dat_rx(i1), dat_ry(i1), 'o');
    %                 text(dat_rx(i1), dat_ry(i1), ...
    %                      sprintf('%s:%s', repProjs{i1, 1}, repProjs{i1, 2}), ...
    %                      'Color', 'b', 'FontSize', 7);
                end            
    
                xs = get(gca, 'XLim');
                ys = get(gca, 'YLim');
                lims = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
                set(gca, 'XLim', lims, 'YLim', lims);
                grid on;

                xlabel(sprintf('%s (Session 1)', measName));
                ylabel(sprintf('%s (Session 2)', measName));
                axis square;

                plot(lims, lims, '--', 'Color', [0.5, 0.5, 0.5]);


                xs = get(gca, 'XLim');
                ys = get(gca, 'YLim');

                text(xs(1) + 0.05 * (xs(2) - xs(1)), ...
                     ys(2) - 0.06 * (ys(2) - ys(1)), ...
                     sprintf('lincorr: R^2=%f; p=%f', lc_r2, lc_p));
                plot(xs, xs * lc_k(2) + lc_k(1), 'b-');


                text(xs(1) + 0.05 * (xs(2) - xs(1)), ...
                     ys(2) - 0.12 * (ys(2) - ys(1)), ...
                     sprintf('spear: rho=%f; p=%f', sp_r, sp_p));

                set(gca, 'XLim', xs, 'YLim', ys);
            end
        elseif length(size(dat)) == 3
            for i1 = 1 : size(repTab, 1)
                figure('Name', sprintf('SubjID: %s / %s', sids{repTab(i1, 1)}, sids{repTab(i1, 2)}), ...
                       'Position', [100, 200, 900, 900], 'Color', [1, 1, 1]);
                for i2 = 1 : 2
                    if i2 == 1
                        dat_rx = dat(:, :, repTab(i1, 1));
                        dat_ry = dat(:, :, repTab(i1, 2));
                        clr = 'b';

                        mats{1} = dat_rx;
                        mats{2} = dat_ry;

                    else
                        dat_rx = dat(:, :, repTab(i1, 1));
                        idxOther = setxor(1 : size(dat, 3), repTab(i1, 1));
                        idxOther = setxor(idxOther, repTab(i1, 2));

                        if args.bRandOther
                            rpOthers = randperm(length(idxOther));
                            sidComp = sids{idxOther(rpOthers(1))};
                            dat_ry = mean(dat(:, :, idxOther(rpOthers(1))), 3);
                        else
                            dat_ry = mean(dat(:, :, idxOther), 3);
                        end

                        clr = 'r';                    
                        mats{3} = dat_ry;
                    end

                    dat_rx = dat_rx(:);
                    dat_ry = dat_ry(:);

                    subplot(2, 2, i2);
                    hold on;
                    plot(dat_rx, dat_ry, 'o', 'Color', clr);

                    set(gca, 'XScale', 'log', 'YScale', 'log');
                    axis square;

                    xs = get(gca, 'XLim');
                    ys = get(gca, 'YLim');
                    lims = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
                    if lims(1) == 0
                        lims(1) = lims(2) / 1e2;
                    end
                    set(gca, 'XLim', lims, 'YLim', lims);
                    grid on;
                    plot(lims, lims, '-', 'Color', [0.5, 0.5, 0.5]);

                    xlabel(sprintf('%s from Session 1: %s (%s)', strrep(args.anaType, '_', '\_'), ...
                           strrep(sids{repTab(i1, 1)}, '_', '\_'), projs{repTab(i1, 1)}));
                    if i2 == 1
                        ylabel(sprintf('%s from Session 2: %s (%s)', strrep(args.anaType, '_', '\_'), ...
                               strrep(sids{repTab(i1, 2)}, '_', '\_'), projs{repTab(i1, 2)}));
                        titl = 'Session 2 vs. session 1';
                    else
                        if args.bRandOther
                            ylabel(sprintf('Session from another subject: %s', sidComp));
                            titl = sprintf('Another subject (%s) vs. session 1', sidComp);
                        else
                            ylabel(sprintf('%s from average of %d other subjects', ...
                                   strrep(args.anaType, '_', '\_'), length(idxOther)));
                            titl = 'Other-average vs. session 1';
                        end
                    end

                    %-- calculate correlation coeffieicents and differences --%
                    dat_rx(dat_rx == 0) = MIN_TNS;
                    dat_ry(dat_ry == 0) = MIN_TNS;

                    [rho_spear, t_spear, p_spear] = spear(log(dat_rx), log(dat_ry));
                    text(xs(1) + 1e-8 * range(xs), ys(2) - 0.06 * range(ys), ...
                         sprintf('rho = %f, p = %e', rho_spear, p_spear));
                    title(titl);

                    [k_lin, r2_lin, p_lin] = lincorr(log(dat_rx), log(dat_ry));
                    text(xs(1) + 1e-8 * range(xs), ys(2) - 0.5 * range(ys), ...
                         sprintf('r^2 = %f, p = %e', r2_lin, p_lin));

                    mean_diff = mean(abs(log(dat_rx) - log(dat_ry)));
                    sd_diff = std(abs(log(dat_rx) - log(dat_ry)));
                    text(xs(1) + 1e-8 * range(xs), ys(2) - 0.9 * range(ys), ...
                         sprintf('Mean diff. = %f; SD diff. = %f', mean_diff, sd_diff));
                end

                for i2 = 1 : 3
                    subplot('Position', [0.05 + (i2 - 1) * 0.3, 0.15, 0.28, 0.325]);
                    imagesc(mats{i2});

                    if i2 == 1
                        title(sprintf('Data from Session 1: %s (%s)', ...
                              strrep(sids{repTab(i1, 1)}, '_', '\_'), projs{repTab(i1, 1)}));
                    elseif i2 == 2
                        title(sprintf('Data from Session 2: %s (%s)', ...
                              strrep(sids{repTab(i1, 2)}, '_', '\_'), projs{repTab(i1, 2)}));
                    else
                        if args.bRandOther
                            title(sprintf('Data from another subject: %s', sidComp));                       
                        else
                            title(sprintf('Average from average of %d subjects', length(sids)));                        
                        end
                    end
                end
            end
        
        end
    end
    
elseif isequal(args.grpScheme, 'byGroup') %-- Between group comparisons --%
    %-- Comparison, project by project --%
    ng = max(grp);
    
    figure('Position', [100, 100, 1200, 600], 'Color', 'w');
    for i1 = 1 : numel(uprojs)
        subplot(2, 3, i1);
        hold on;
        title(uprojs{i1});
        
        t_grp = grp(fsic(projs, uprojs{i1}));        
        t_dat = dat(fsic(projs, uprojs{i1}));        
        
        c_dat = {};
        for i2 = 1 : ng
            gdat = t_dat(t_grp == i2);
            if ~isempty(gdat)
                c_dat{end + 1} = gdat;
            end
            
            plot(i2, nanmean(gdat), 'o');
            plot([i2, i2], nanmean(gdat) + [-1, 1] * nanste(gdat), '-');
            plot(repmat(i2 + 0.2, 1, length(gdat)), gdat, 'o');
        end
        
        set(gca, 'XLim', [0, ng + 1], 'XTick', 1 : ng, 'XTIckLabel', args.groups);
        
        xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
        strN = 'N = (';
        for i2 = 1 : numel(c_dat)
            strN = sprintf('%s%d, ', strN, length(c_dat{i2}));            
        end
        strN = [strN(1 : end - 2), ')'];
        text(xs(1) + 0.05 * range(xs), ys(2) - 0.05 * range(ys), strN);
        
        if length(c_dat) == 2 && ...
           length(c_dat{1}) > 2 && length(c_dat{2}) > 2 
           %-- Perform t-test and ranksum test --%
            [t_h, t_p, t_ci, t_stats] = ttest2(c_dat{1}, c_dat{2});
            
            rs_p = ranksum(c_dat{1}, c_dat{2});
            
            text(xs(1) + 0.05 * range(xs), ys(2) - 0.10 * range(ys), ...
                 sprintf('t-test: t=%f, p=%f', t_stats.tstat, t_p));
            text(xs(1) + 0.05 * range(xs), ys(2) - 0.15 * range(ys), ...
                 sprintf('ranksum: p=%f', rs_p))
        elseif length(c_dat) > 2 % Perform ANOVA
            
        end
    end
    
    %-- Perform two-way ANOVA, generally unbalenced --%
    if length(uprojs) > 1
        aovTab = anovan(dat, {grp, projs}, 'display', 'off');
        info_log('');
        info_log(sprintf('Two-way ANOVA: '))
        info_log(sprintf('\tGROUP: p = %f', aovTab(1)))
        info_log(sprintf('\tPROJ: p = %f', aovTab(2)));
    end
    
end

return