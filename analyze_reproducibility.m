function analyze_reproducibility()
%% CONFIG
WM_DEPTH = 2;
HEMIS = {'lh', 'rh'};

%% Get the list of all speech ROIs
parc = load('aparc12.mat');
assert(exist('parc', 'var') == 1);
assert(isfield(parc, 'speechROIs'));

rois = cell(1, size(parc.speechROIs, 1));
for i1 = 1 : size(parc.speechROIs, 1)
    rois{i1} = deblank(parc.speechROIs(i1, :));
end

%% TNS reproducibility 
nROIs = length(rois);
nTracts = nROIs * (nROIs - 1) / 2;
res.TNS.lc_r = nan(1, nTracts * 2);
res.TNS.sp_r = nan(1, nTracts * 2);

rc = 1;
for i1 = 1 : length(HEMIS)
    hemi = HEMIS{i1};
    
    for i2 = 1 : length(rois)
        for i3 = i2 + 1 : length(rois)
            tract = sprintf('%s_%s-%s_%s', hemi, rois{i2}, hemi, rois{i3});
            t_res = dwi_group('TNS', 'rep', '--roi', tract, '--symm-tract');
            
            res.TNS.lc_r(rc) = t_res(1);
            res.TNS.sp_r(rc) = t_res(2);
            rc = rc + 1;
        end
    end
end

%% FA reproducibility
res = struct;
res.FA.lc_r = nan(1, length(rois) * 2);
res.FA.sp_r = nan(1, length(rois) * 2);

for i1 = 1 : length(HEMIS)
    hemi =  HEMIS{i1};
    for j1 = 1 : length(rois)
        t_roi = strcat(hemi, '_', rois{j1});
        t_res = dwi_group('FA', 'rep', '--roi', t_roi, '--wm-depth', WM_DEPTH);
        
        res.FA.lc_r((i1 - 1) * length(rois) + j1) = t_res(1);
        res.FA.sp_r((i1 - 1) * length(rois) + j1) = t_res(2);
    end
end



return