function dat = get_tensor_measure(tmFile, meas, roi, wmDepth, varargin)
%%
% roi == all: get the arithmetic mean of all regions in the mat file
% Option --norm-by-all: normize the value by the mean of all regions

%%
tmdat = load(tmFile);

%% Optional input arguments
bNormByAll = ~isempty(fsic(varargin, '--norm-by-all'));

%% Load the depth
jd = find(tmdat.wmDepths == wmDepth);
if length(jd) ~= 1
    error_log(sprintf('Cannot find WM depth %d in tensor measure file: %s', ...
                      wmDepth, tmFile));
end

%% Locate ROI
if ~isequal(roi, 'all')
    bROIFound = 0;
    for i1 = 1 : size(tmdat.roiNames, 1)
        t_roi = deblank(tmdat.roiNames(i1, :));
    
        if isequal(t_roi, roi)
            bROIFound = 1;
            break;
        end
    end

    if ~bROIFound
        error_log(sprintf('Cannot find ROI %s in tensor measure file: %s', ...
                          roi, tmFile));
    end
    jr = i1;
end
    
%%
meanAll = mean(tmdat.tensMeas.(meas)(:, jd));
if ~isequal(roi, 'all')
    dat = tmdat.tensMeas.(meas)(jr, jd);
    if bNormByAll
        dat = dat / meanAll;
    end
else
    dat = meanAll;
end

return