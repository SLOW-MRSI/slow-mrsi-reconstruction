function plotBrainMaskPreview(brainMask, previewMode, waterMapSmooth, params)
%PLOTBRAINMASKPREVIEW Optional quick preview of the generated brain mask.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

if ~strcmpi(previewMode, 'on')
    return
end

figure(4); clf
set(gcf, 'Name', 'Brain mask adjustment', 'NumberTitle', 'off', ...
    'ToolBar', 'figure', 'MenuBar', 'figure')
hasWaterMap = nargin >= 3 && ~isempty(waterMapSmooth);
nSlices = size(brainMask, 3);
nRows = 2;
nCols = ceil(nSlices / nRows);
left = 0.035;
right = 0.985;
bottom = 0.13;
top = 0.875;
horizontalGap = 0.028;
verticalGap = 0.08;
axisWidth = (right - left - (nCols - 1) * horizontalGap) / nCols;
axisHeight = (top - bottom - (nRows - 1) * verticalGap) / nRows;

for nz = 1:nSlices
    row = floor((nz - 1) / nCols) + 1;
    col = mod(nz - 1, nCols) + 1;
    axisLeft = left + (col - 1) * (axisWidth + horizontalGap);
    axisBottom = top - row * axisHeight - (row - 1) * verticalGap;
    axes('Position', [axisLeft axisBottom axisWidth axisHeight])
    if hasWaterMap
        imagesc(waterMapSmooth(:,:,nz))
        axis image off
        colormap(gca, gray)
        hold on
        if any(brainMask(:,:,nz), 'all')
            contour(brainMask(:,:,nz), [0.5 0.5], 'r', 'LineWidth', 0.8)
        end
        hold off
    else
        imshow(brainMask(:,:,nz))
    end
    title(sprintf('z = %d', nz), 'FontSize', 9, 'FontWeight', 'bold')
end

if nargin >= 4 && isstruct(params)
    annotation('textbox', [0.02 0.905 0.96 0.04], ...
        'String', sprintf('Brain mask: threshold %.4g, skull margin %d voxels', ...
        params.thresholdFraction, params.skullMarginVoxels), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'FontSize', 11, 'FontWeight', 'bold');
    annotation('textbox', [0.02 0.955 0.96 0.025], ...
        'String', 'Inspect the water map and red mask boundary, then update parameters or accept.', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'FontSize', 8);
end
end
