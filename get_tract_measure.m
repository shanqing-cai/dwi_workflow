function dat = get_tract_measure(connFile, anaType, ...
                                 seedHemi, seedROI, targHemi, targROI)
%%

%%
cdat = load(connFile);

%% Locate the ROIs
idxROIs = nan(1, 2);

assert(~isempty(seedHemi));
assert(~isempty(seedROI));

if isempty(targHemi) && isempty(targROI)
    nROIs = 1;
else
    nROIs = 2;
end

for i1 = 1 : nROIs
    if i1 == 1
        roi = seedROI;
    else
        roi = targROI;
    end
    
    bFound = 0;
    for i2 = 1 : size(cdat.roiNames, 1)
        t_roi = deblank(cdat.roiNames(i2, :));
        
        if isequal(t_roi, roi)
            bFound = 1;
            break;
        end        
    end
    
    if ~bFound
        error_log(sprintf('Cannot find ROI %s in connectivity mat file: %s', ...
                          roi, connFile));
    end
    idxROIs(i1) = i2;    
end

%% Get data
if nROIs == 2
    dat = cdat.connMat(idxROIs(1), idxROIs(2));
else
    dat = cdat.connMat(idxROIs(1), :);
end
return