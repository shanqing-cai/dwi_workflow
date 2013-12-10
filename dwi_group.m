function dwi_group(anaType, grpScheme, varargin)
%% dwi_group
% Input arguments:
%       anaType: analysis type {FA, MD, TNS}
%                TNS - tract density normalized by seed size
%       grpScheme: grouping scheme {byStudy, byGroup, rep, byAge}
%                   rep - Reproducibility test: multiple sessions on the
%                         same subject
%                   age - Analayze age effects
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
%       --no-self:     remove the self-projections in the TNS analysis
%                      Self-projections refers to the projection from an ROI to itself                    
%       --rand-other: Use random other individual subject to compare with
%                     the TNS measure of the repeated subject.
%       -v | --verbose: verbose mode
%
%%
ANALYSIS_TYPES = {'FA', 'MD', 'TNS'};
GROUPING_SCHEMES = {'byStudy', 'byGroup', 'rep', 'byAge'};

analysisSettingsMat = 'dwi_analysis_settings.mat';
projInfoMat = 'dwi_project_info.mat';

DEFAULT_PARC = 'aparc12';

%% Other constants
MIN_TNS = 1e-6;

%--- Visualization settings ---%
fontSize = 15;

%% Visualization options
COLORS = {[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0.5, 1], ...
          [0.5, 0.25, 0], [0.5, 0, 1], [0.5, 0.5, 0], [0, 1, 1], ...
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

args.bNormByAll = ~isempty(fsic(varargin, '--norm-by-all'));
args.bNoSelf = ~isempty(fsic(varargin, '--no-self'));

args.bRandOther = ~isempty(fsic(varargin, '--rand-other'));

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
            else
                [t_dat, t_roiNames] = get_tract_measure(connFile, args.anaType, ...
                                                        seedHemi, 'all');
                if ~isempty(prev_roiNames)
                    assert(isequal(t_roiNames, prev_roiNames));
                end
            end
        else
            error_log(sprintf('Analaysis type %s has not been implemented yet', ...
                              args.anaType));
        end
        
        if isequal(args.roi, 'all') || isequal(args.roi, 'lh_all') || isequal(args.roi, 'rh_all')
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

    
%% Some data formatting

%% Statistical analysis and visualization
measName = strrep(sprintf('%s: %s', args.anaType, args.roi), '_', '\_');
if args.bNormByAll
    measName = [measName, ' (norm. by allAvg)'];
end

if isequal(args.grpScheme, 'byStudy') || isequal(args.grpScheme, 'byAge')
    ugrps = unique(grp);

    fhdl = figure('Name', sprintf('%s: %s', args.grpScheme, args.anaType), 'Color', 'w');
    set(gca, 'FontSize', fontSize);
    hold on; box on;
    
    if length(size(dat)) == 2 && isvector(dat) %--- Scalar measure from each subject ---%
        [aov1_p, aov1_table] = anova1(dat, grp, 'off');
        for i1 = 1 : numel(ugrps)
            t_grp = ugrps(i1);
            idx = find(grp == t_grp);

            if length(idx) > 0
                if isequal(args.grpScheme, 'byStudy')
                    plot(idx, dat(idx), 'o-', 'Color', COLORS{i1});
                    text(idx(1) + 0.5, max(dat(idx)), ...
                         deblank(projInfo.name(t_grp, :)), ...
                         'FontSize', fontSize, 'Color', COLORS{i1});
                elseif isequal(args.grpScheme, 'byAge')
                    plot(ages(idx), dat(idx), 'o', 'Color', COLORS{i1});
                    text(ages(idx(1)), dat(idx(1)), ...
                         deblank(projInfo.name(t_grp, :)), ...
                         'FontSize', fontSize, 'Color', COLORS{i1});
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
    if length(size(dat)) >= 2
        legend(t_legend);
        
    end
   
elseif isequal(args.grpScheme, 'rep')
    if isequal(args.anaType, 'FA') || isequal(args.anaType, 'MD')
        if isvector(dat) % Single ROI
            
        else % A set of ROIs;
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
        end
        
        figure('Color', 'w'); 
        set(gca, 'FontSize', fontSize)
        hold on;
        plot(r2s_self, r2s_so, 'ko');
        lims = [0, 1];
        set(gca, 'XLim', lims, 'YLim', lims);
        grid on;
        plot(lims, lims, '-', 'Color', [0.5, 0.5, 0.5]);
        xlabel('R^2 of self-self correlation');
        ylabel('R^2 of self-other correlation');       
        
    elseif isequal(args.anaType, 'TNS')
    
        if length(size(dat)) == 2   % Single connection
            dat_rx = dat(repTab(:, 1));
            dat_ry = dat(repTab(:, 2));

            figure;
            hold on;
            for i1 = 1 : length(dat_rx) 
                plot(dat_rx(i1), dat_ry(i1), 'o');
                text(dat_rx(i1), dat_ry(i1), ...
                     sprintf('%s:%s', repProjs{i1, 1}, repProjs{i1, 2}), ...
                     'Color', 'b', 'FontSize', 7);
            end

            xs = get(gca, 'XLim');
            ys = get(gca, 'YLim');
            lims = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
            set(gca, 'XLim', lims, 'YLim', lims);
            grid on;

            xlabel(sprintf('%s (Session 1)', measName));
            ylabel(sprintf('%s (Session 2)', measName));
            axis square;
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
end

return