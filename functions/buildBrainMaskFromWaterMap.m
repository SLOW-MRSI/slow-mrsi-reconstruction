function [brainMask, mask_main, mask_wat, thresholdValue, skullMarginVoxels] = buildBrainMaskFromWaterMap( ...
    waterMapSmooth, thresholdFraction, skullMarginVoxels, minObjectVoxels, closeRadius)
%BUILDBRAINMASKFROMWATERMAP Build an inner brain mask from smoothed water signal.
% The mask is made from a filled slice-wise head outline, then shrunk inward
% by a manually selected skull margin.

if nargin < 3 || isempty(skullMarginVoxels)
    skullMarginVoxels = 3;
end
if nargin < 4 || isempty(minObjectVoxels)
    minObjectVoxels = 100;
end
if nargin < 5 || isempty(closeRadius)
    closeRadius = 2;
end

thresholdFraction = max(0, double(thresholdFraction));
skullMarginVoxels = max(0, round(double(skullMarginVoxels)));
closeRadius = max(0, round(double(closeRadius)));

maxWaterSignal = max(waterMapSmooth(:));
if isempty(maxWaterSignal) || ~isfinite(maxWaterSignal) || maxWaterSignal <= 0
    thresholdValue = 0;
    mask_wat = false(size(waterMapSmooth));
    mask_main = mask_wat;
    brainMask = mask_wat;
    return
end

thresholdValue = thresholdFraction * maxWaterSignal;
mask_wat = waterMapSmooth > thresholdValue;

mask_main = false(size(mask_wat));
minSliceVoxels = max(5, round(minObjectVoxels / max(1, size(mask_wat, 3))));
for kz = 1:size(mask_wat, 3)
    sliceMask = mask_wat(:,:,kz);
    if closeRadius > 0
        sliceMask = imclose(sliceMask, strel('disk', closeRadius, 0));
    end
    sliceMask = imfill(sliceMask, 'holes');
    sliceMask = bwareaopen(sliceMask, minSliceVoxels, 8);

    CC = bwconncomp(sliceMask, 8);
    if CC.NumObjects > 0
        numPix = cellfun(@numel, CC.PixelIdxList);
        [~, idLargest] = max(numPix);
        sliceMain = false(size(sliceMask));
        sliceMain(CC.PixelIdxList{idLargest}) = true;
        mask_main(:,:,kz) = sliceMain;
    end
end

brainMask = false(size(mask_main));
for kz = 1:size(mask_main, 3)
    if skullMarginVoxels > 0
        distanceFromBoundary = bwdist(~mask_main(:,:,kz));
        brainMask(:,:,kz) = distanceFromBoundary > skullMarginVoxels;
    else
        brainMask(:,:,kz) = mask_main(:,:,kz);
    end
end
end
